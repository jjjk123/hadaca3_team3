##################################################################################################
### PLEASE only edit the program function between YOUR CODE BEGINS/ENDS HERE                   ###
##################################################################################################


########################################################
### Package dependencies /!\ DO NOT CHANGE THIS PART ###
########################################################
import subprocess
import sys
import importlib

def program(mix_rna=None, ref_rna=None, mix_met=None, ref_met=None, **kwargs):

  ##
  ## YOUR CODE BEGINS HERE
  ##

  required_packages = ["sklearn","pandas",'scipy']
  install_and_import_packages(required_packages)
  from sklearn.linear_model import LinearRegression
  from sklearn.preprocessing import StandardScaler
  import pandas
  import scipy
  import numpy

  from attachement.additionnal_script import find_optimal_weights, calculate_weighted_sum
  #additionnal_script.useless_function()

  meth_rna_mapping = pandas.read_csv("attachement/mapping_meth_rna.csv", index_col=0)
  common_bulk_meth = meth_rna_mapping["UCSC_RefGene_Name"].values

  filtered_mix_rna = mix_rna.loc[common_bulk_meth]
  filtered_mix_met = mix_met.loc[meth_rna_mapping.index]
  filtered_mix_met.index = meth_rna_mapping["UCSC_RefGene_Name"].values
  filtered_ref_rna = ref_rna.loc[common_bulk_meth]
  filtered_ref_meth = ref_met.loc[meth_rna_mapping.index]
  filtered_ref_meth.index = meth_rna_mapping["UCSC_RefGene_Name"].values

  pseudo_bulk_peng = pandas.read_csv("attachement/peng_pseudo_bulk_sum.csv", index_col=0)
  pseudo_bulk_peng_cpm = pseudo_bulk_peng/pseudo_bulk_peng.sum(axis=0)*1000000
  pseudo_bulk_peng_cpm = pseudo_bulk_peng_cpm.loc[pseudo_bulk_peng_cpm.sum(axis=1)>0]

  global_common_genes = numpy.intersect1d(pseudo_bulk_peng_cpm.index, common_bulk_meth)

  def estimate_proportions(mix_df, ref_df):
    results = []
    for i in range(len(mix_df.columns)):
        mix_col = mix_df.iloc[:, i]  # Select the i-th column as a Series
        #res = LinearRegression(fit_intercept=False).fit(ref_df, mix_col).coef_
        res, _ = scipy.optimize.nnls(ref_df.to_numpy(), mix_col.to_numpy())
        res[res < 0] = 0
        results.append(res)

    # Normalize the results to get proportions
    props = pandas.DataFrame([res_i / sum(res_i) for res_i in results], columns=ref_df.columns)
    return props.T
  
  # Filter CpG sites
  row_max = numpy.max(filtered_ref_meth.values, axis=1)
  row_means_without_max = (numpy.sum(filtered_ref_meth.values, axis=1) - row_max) / (filtered_ref_meth.values.shape[1] - 1)
  diffs = row_max - row_means_without_max
  threshold = diffs[diffs>numpy.min(numpy.mean(filtered_ref_meth.values, axis=0))]
  sorted_indices = numpy.argsort(-diffs)
  sorted_indices_reduced = sorted_indices[:threshold.shape[0]]
  #filtered_ref_meth = filtered_ref_meth.iloc[sorted_indices_reduced]


  # Creation of an index, idx_feat, corresponding to the intersection of features present in the references and those present in the mixtures.
  idx_feat_meth = filtered_mix_met.index.intersection(filtered_ref_meth.index)
  mix_filtered_meth = filtered_mix_met.loc[idx_feat_meth, :]
  ref_filtered_meth = filtered_ref_meth.loc[idx_feat_meth, :]

  # Scaling does not work
  prop_meth = estimate_proportions(mix_filtered_meth, ref_filtered_meth)
  # Labeling of estimated proportions 
  prop_meth.columns = mix_filtered_meth.columns

  # Creation of an index, idx_feat, corresponding to the intersection of features present in the references and those present in the mixtures.
  #idx_feat_rna = filtered_mix_rna.index.intersection(filtered_ref_rna.index)
  mix_filtered_rna = filtered_mix_rna.loc[global_common_genes, :]
  ref_filtered_rna = filtered_ref_rna.loc[global_common_genes, :]

  prop_rna = estimate_proportions(mix_filtered_rna, ref_filtered_rna)

  # Labeling of estimated proportions 
  prop_rna.columns = mix_filtered_rna.columns

  # Creation of an index, idx_feat, corresponding to the intersection of features present in the references and those present in the mixtures.
  #idx_feat_pseudo_rna = filtered_mix_rna.index.intersection(pseudo_bulk_peng_cpm.index)
  mix_filtered_pseudo_rna = filtered_mix_rna.loc[global_common_genes, :]
  ref_filtered_pseudo_rna = pseudo_bulk_peng_cpm.loc[global_common_genes, :]

  prop_pseudo_rna = estimate_proportions(mix_filtered_pseudo_rna, ref_filtered_pseudo_rna)

  # Labeling of estimated proportions 
  prop_pseudo_rna.columns = mix_filtered_rna.columns

  # Just summing is not the best
  #final_prop = (prop_rna + prop_meth)/2

  # Multiplication of RNA and Meth doesn't work
#   mix_rna_times_meth = (1-filtered_mix_met).mul(filtered_mix_rna)
#   ref_rna_times_meth = (1-ref_filtered_meth).mul(ref_filtered_rna)
#   prop_rna_times_meth = estimate_proportions(mix_rna_times_meth, ref_rna_times_meth)
#   prop_rna_times_meth.columns = mix_rna_times_meth.columns

  y_hat = [ref_filtered_pseudo_rna.dot(prop_pseudo_rna), ref_filtered_rna.dot(prop_rna), ref_filtered_meth.dot(prop_meth)]
  best_weights, best_rmse = find_optimal_weights((mix_filtered_rna+mix_filtered_meth)/2, y_hat)
  print(best_weights, best_rmse)

  p_hat = [prop_pseudo_rna, prop_rna, prop_meth]
  final_prop = pandas.DataFrame(calculate_weighted_sum(p_hat, best_weights),
                                index=prop_pseudo_rna.index, columns=prop_pseudo_rna.columns)

  return final_prop
  ##
  ## YOUR CODE ENDS HERE
  ##


# Install and import each package
def install_and_import_packages(required_packages):
  for package in required_packages:
      try:
          globals()[package] = importlib.import_module(package)
      except ImportError:
          print('impossible to import, installing packages',package)
          package_to_install = 'scikit-learn' if package == 'sklearn' else package
          subprocess.check_call([sys.executable, "-m", "pip", "install", package_to_install])
          globals()[package] = importlib.import_module(package)

def validate_pred(pred, nb_samples=None, nb_cells=None, col_names=None):
    error_status = 0  # 0 means no errors, 1 means "Fatal errors" and 2 means "Warning"
    error_informations = ''

    # Ensure that all sum of cells proportion approximately equal 1
    if not numpy.allclose(numpy.sum(pred, axis=0), 1):
        msg = "The prediction matrix does not respect the laws of proportions: the sum of each column should be equal to 1\n"
        error_informations += msg
        error_status = 2

    # Ensure that the prediction has the correct names
    if not set(col_names) == set(pred.index):
        msg = f"The row names in the prediction matrix should match: {col_names}\n"
        error_informations += msg
        error_status = 2

    # Ensure that the prediction returns the correct number of samples and number of cells
    if pred.shape != (nb_cells, nb_samples):
        msg = f'The prediction matrix has the dimension: {pred.shape} whereas the dimension: {(nb_cells, nb_samples)} is expected\n'
        error_informations += msg
        error_status = 1

    if error_status == 1:
        # The error is blocking and should therefore stop the execution
        raise ValueError(error_informations)
    if error_status == 2:
        print("Warning:")
        print(error_informations)


##############################################################
### Generate a prediction file /!\ DO NOT CHANGE THIS PART ###
##############################################################

# List of required packages
required_packages = [
  "numpy",
  "pandas",
  "rpy2",
  "zipfile",
  "inspect",
]
install_and_import_packages(required_packages)

# from rpy2.robjects import pandas2ri
import os
import rpy2.robjects
readRDS = rpy2.robjects.r['readRDS']
saveRDS= rpy2.robjects.r["saveRDS"]

from rpy2.robjects import pandas2ri
pandas2ri.activate()


r_code_get_rowandcolnames = '''
get_both <- function(ref_names = "reference_fruits.rds", mat = NULL) {
  ref_names <- readRDS(ref_names)
  if (!is.null(mat)) {
    return(list(  rownames(ref_names[[mat]]), colnames(ref_names[[mat]]) ))
  } else {
    return(list(rownames(ref_names),colnames(ref_names) ))
  }
}
'''
rpy2.robjects.r(r_code_get_rowandcolnames)
get_both_row_col = rpy2.robjects.r['get_both']


# # Function to convert R object to pandas DataFrame or numpy array
def r_object_to_python(r_object,file,element_name):
    try:
        # Try to convert to pandas DataFrame
        return pandas2ri.rpy2py(r_object)
    except NotImplementedError:
        rows, columns =get_both_row_col(file,element_name)
        if(isinstance(columns, type (rpy2.robjects.NULL))):
            df = pandas.DataFrame(r_object, index=rows)
        else: 
            columns = list(columns)
            df = pandas.DataFrame(r_object, columns=columns, index=rows)
        return df

# Function to extract named data elements and convert appropriately
def extract_data_element(data, file, element_name):
    if element_name in data.names:
        element = data.rx2(element_name)
        return r_object_to_python(element,file,element_name)
    return None


dir_name = "data"+os.sep

datasets_list = [filename for filename in os.listdir(dir_name) if filename.startswith("mixes")]

ref_file = os.path.join(dir_name, "reference_pdac.rds")
print("reading reference file")
reference_data = readRDS(ref_file)
ref_bulkRNA = extract_data_element(reference_data,ref_file, 'ref_bulkRNA') 
ref_met = extract_data_element(reference_data,ref_file, 'ref_met')
ref_scRNA = extract_data_element(reference_data,ref_file, 'ref_scRNA')

predi_dic = {}
for dataset_name in datasets_list:

    file= os.path.join(dir_name,dataset_name)
    mixes_data = readRDS(file)

    print(f"generating prediction for dataset: {dataset_name}")

    mix_rna = extract_data_element(mixes_data,file, 'mix_rna') 
    mix_met = extract_data_element(mixes_data,file, 'mix_met')

    pred_prop = program(mix_rna, ref_bulkRNA, mix_met=mix_met, ref_met=ref_met)
    validate_pred(pred_prop, nb_samples=mix_rna.shape[1], nb_cells=ref_bulkRNA.shape[1], col_names=ref_bulkRNA.columns)
    predi_dic[dataset_name] = pred_prop

############################### 
### Code submission mode

# we generate a zip file with the 'program' source code

if not os.path.exists("submissions"):
    os.makedirs("submissions")

# we save the source code as a Python file named 'program.py':
with open(os.path.join("submissions", "program.py"), 'w') as f:
    f.write(inspect.getsource(program))

date_suffix = pandas.Timestamp.now().strftime("%Y_%m_%d_%H_%M_%S")




# we create the associated zip file:
zip_program = os.path.join("submissions", f"program_{date_suffix}.zip")
with zipfile.ZipFile(zip_program, 'w') as zipf:
    zipf.write(os.path.join("submissions", "program.py"), arcname="program.py")


def zipdir(path, ziph):
    # ziph is zipfile handle
    for root, dirs, files in os.walk(path):
        for file in files:
            ziph.write(os.path.join(root, file), 
                       os.path.relpath(os.path.join(root, file), 
                                       os.path.join(path, '..')))
if os.path.exists("attachement"):
    with zipfile.ZipFile(zip_program, 'a', zipfile.ZIP_DEFLATED) as zipf:
        zipdir('attachement/', zipf)



# # Check if the "attachment" directory exists
# if os.path.exists("attachement"):
#     # Append the contents of the "attachment" directory to the zip archive
#     with zipfile.ZipFile(zip_program, mode="a") as zf:
#         for root, _, files in os.walk("attachement"):
#             for file in files:
#                 file_path = os.path.join(root, file)
#                 # Add file to zip while preserving directory structure
#                 arcname = os.path.relpath(file_path, start="attachement")
#                 zf.write(file_path, arcname)

print(zip_program)

###############################
### Result submission mode  

# Generate a zip file with the prediction
if not os.path.exists("submissions"):
    os.makedirs("submissions")

prediction_name = "prediction.rds"

saveRDS(rpy2.robjects.ListVector(predi_dic), os.path.join("submissions", prediction_name))



# Create the associated zip file:
zip_results = os.path.join("submissions", f"results_{date_suffix}.zip")
with zipfile.ZipFile(zip_results, 'w') as zipf:
    zipf.write(os.path.join("submissions", prediction_name), arcname=prediction_name)

print(zip_results)

