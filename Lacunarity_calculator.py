import time
import sys

def type_text(text, delay=0.1):
    for char in text:
        sys.stdout.write(char)
        sys.stdout.flush()
        time.sleep(delay)

type_text("Initiating Code....")

try:
    import imea
    from imea.measure_2d.micro import fractal_dimension_boxcounting as fract_dim
    from scipy.signal import convolve2d as conv2d
    import multiprocessing as mp    
    import cv2
    import nibabel as nb 
    import pandas as pd
    import numpy as np   
    import os      
    import warnings
    warnings.filterwarnings('ignore')
except Exception as e:
    print(f"An error occurred while importing a package: {e}. Check the package list that needs to be installed.")



t0 = time.time()


#Lacunarity Calculation function - Divided into 3 functions for each subcomponent.
def lac1(x):                     #non enhancing lac   
   bin_ = np.zeros((x.shape[0], x.shape[1]))   
   bin_[x!=0] = 1   
   bin_ = bin_.astype(np.uint8)
   
   contours = cv2.findContours(image = bin_, mode = cv2.RETR_TREE, method = cv2.CHAIN_APPROX_SIMPLE)[0]
   contours = sorted(contours, key = cv2.contourArea, reverse = True)
   
   if len(contours) == 0:
      return 0
   
   c_0 = contours[0]              # Get the 4 points of the bounding rectangle
   xx, yy, w, h = cv2.boundingRect(c_0)
   X  =  x[yy:yy+h, xx:xx+w]
   box_size =  [i for i in range(2, 240) if ((i)**2 < (0.45)*(X.shape[0]* X.shape[1])) and i <= min(X.shape[0], X.shape[1])]
  
   LAMBDA = []  
   
   XX = np.zeros((X.shape[0], X.shape[1]))   
   XX[X==1] = 1

   for box in box_size:
      count = []
      count, edge= np.histogram(np.ravel(conv2d(XX, np.ones((box, box)), mode = 'valid')), bins = [i for i in range(0, (box**2) + 2)])#(box**2) + 1)
      q = count/(XX.shape[0] - box + 1)**2
      s = np.array([i for i in range(0,box**2 + 1)])
      lam_bda =  sum((s**2)*q)/(sum(q*s))**2   
      LAMBDA.append(lam_bda)

        
   return np.nanmean(LAMBDA)


def lac2(x):   
   bin_ = np.zeros((x.shape[0], x.shape[1]))  
   bin_[x!=0] = 1   
   bin_ = bin_.astype(np.uint8)  

   contours = cv2.findContours(image = bin_, mode = cv2.RETR_TREE, method = cv2.CHAIN_APPROX_SIMPLE)[0]
   contours = sorted(contours, key = cv2.contourArea, reverse = True)
   
   if len(contours) == 0:
      return 0
   
   c_0 = contours[0] # Get the 4 points of the bounding rectangle
   xx, yy, w, h = cv2.boundingRect(c_0)
   X  =  x[yy:yy+h, xx:xx+w]
   
   box_size =  [i for i in range(2, 240) if ((i)**2 < (0.45)*(X.shape[0]* X.shape[1])) and i <= min(X.shape[0], X.shape[1])]  
   LAMBDA = []
   XX = np.zeros((X.shape[0], X.shape[1]))
   
   XX[X==2] = 1

   for box in box_size:
      count = []
      count, edge= np.histogram(np.ravel(conv2d(XX, np.ones((box, box)), mode = 'valid')), bins = [i for i in range(0, (box**2) + 2)])#(box**2) + 1)
      q = count/(XX.shape[0] - box + 1)**2
      s = np.array([i for i in range(0,box**2 + 1)])
      lam_bda =  sum((s**2)*q)/(sum(q*s))**2      
   
      LAMBDA.append(lam_bda)

      
   return np.nanmean(LAMBDA)


def lac4(x):   
   bin_ = np.zeros((x.shape[0], x.shape[1]))  
   bin_[x!=0] = 1   
   bin_ = bin_.astype(np.uint8)  

   contours = cv2.findContours(image = bin_, mode = cv2.RETR_TREE, method = cv2.CHAIN_APPROX_SIMPLE)[0]
   contours = sorted(contours, key = cv2.contourArea, reverse = True)
   
   if len(contours) == 0:
      return 0
   
   c_0 = contours[0]                   # Get the 4 points of the bounding rectangle
   xx, yy, w, h = cv2.boundingRect(c_0)
   X  =  x[yy:yy+h, xx:xx+w]
   
   box_size =  [i for i in range(2, 240) if ((i)**2 < (0.45)*(X.shape[0]* X.shape[1])) and i <= min(X.shape[0], X.shape[1])]
   LAMBDA = []
   XX = np.zeros((X.shape[0], X.shape[1]))
   
   XX[X==4] = 1

   for box in box_size:
      count = []
      count, edge= np.histogram(np.ravel(conv2d(XX, np.ones((box, box)), mode = 'valid')), bins = [i for i in range(0, (box**2) + 2)])#(box**2) + 1)
      q = count/(XX.shape[0] - box + 1)**2
      s = np.array([i for i in range(0,box**2 + 1)])
      lam_bda =  sum((s**2)*q)/(sum(q*s))**2
     
      LAMBDA.append(lam_bda)

      
      
   return np.mean(LAMBDA)

def frac_lac_single_subject(X):   
   X = X.astype(np.uint8)   
   X = np.delete(X, [i for i in range(0, 155) if np.sum(X[:,:,i].ravel()) < 15], 2)
   X_n = np.delete(X, [i for i in range(0, X.shape[-1]) if np.sum(X[:,:,i].ravel() == 1) < 12], 2)
   X_et = np.delete(X, [i for i in range(0, X.shape[-1]) if np.sum(X[:,:,i].ravel() == 4) < 6], 2)
   X_ed = np.delete(X, [i for i in range(0, X.shape[-1]) if np.sum(X[:,:,i].ravel() == 2) < 6], 2)
   
   X_n = [X_n[:,:,i] for i in range(0,X_n.shape[-1]) ]
   X_et = [X_et[:,:,i] for i in range(0,X_et.shape[-1]) ]
   X_ed = [X_ed[:,:,i] for i in range(0,X_ed.shape[-1]) ]
  
   pool1 = mp.Pool(processes = 4)
   pool2 = mp.Pool(processes = 4)
   pool4 = mp.Pool(processes = 4)
     
   nc = np.array(pool1.map(lac1, X_n))
   ed = np.array(pool2.map(lac2, X_ed))
   et = np.array(pool4.map(lac4, X_et))
   
   n_mean =  np.nanmean(nc)
   n_med =  np.nanmedian(nc)
   
   
   ed_mean =  np.nanmean(ed)
   ed_med =  np.nanmedian(ed)
   
   
   et_mean =  np.nanmean(et)
   et_med =  np.nanmedian(et)

   result = np.array([n_mean, n_med, ed_mean, ed_med, et_mean, et_med])
   print('Results: ')
   print('Non Enhancing Lac Mean: ', result[0], '\nEdema Lac Mean: ', result[2], '\nEnhancing Lac Mean: ', result[4])

   return "Code Completed. Exiting....."


def frac_lac(X):   
   X = X.astype(np.uint8)   
   X = np.delete(X, [i for i in range(0, 155) if np.sum(X[:,:,i].ravel()) < 15], 2)
   X_n = np.delete(X, [i for i in range(0, X.shape[-1]) if np.sum(X[:,:,i].ravel() == 1) < 12], 2)
   X_et = np.delete(X, [i for i in range(0, X.shape[-1]) if np.sum(X[:,:,i].ravel() == 4) < 6], 2)
   X_ed = np.delete(X, [i for i in range(0, X.shape[-1]) if np.sum(X[:,:,i].ravel() == 2) < 6], 2)
   
   X_n = [X_n[:,:,i] for i in range(0,X_n.shape[-1]) ]
   X_et = [X_et[:,:,i] for i in range(0,X_et.shape[-1]) ]
   X_ed = [X_ed[:,:,i] for i in range(0,X_ed.shape[-1]) ]
  
   pool1 = mp.Pool(processes = 4)
   pool2 = mp.Pool(processes = 4)
   pool4 = mp.Pool(processes = 4)
     
   nc = np.array(pool1.map(lac1, X_n))
   ed = np.array(pool2.map(lac2, X_ed))
   et = np.array(pool4.map(lac4, X_et))
   
   n_mean =  np.nanmean(nc)
   n_med =  np.nanmedian(nc)
   
   
   ed_mean =  np.nanmean(ed)
   ed_med =  np.nanmedian(ed)
   
   
   et_mean =  np.nanmean(et)
   et_med =  np.nanmedian(et)

   return  np.array([n_mean, n_med, ed_mean, ed_med, et_mean, et_med])

type_text('Done!')
print('')
print('\n*******************')
print('')

print('1. Check code integrity for a single subject.  \n2. Calculate batch lacunarity. \n\nType 1 or 2 to choose.')

choice = str(input())

if choice == '1':
   print('*' * 50)
   print('File path for the image file (Filename ends with: _Glistrboost_Manually_Corrected.nii.gz): ')
   # print('Example: /home/username/directory/file_glistrboost_manually_corrected.nii.gz (or ".nii" file)')
   img_str = str(input())
   print(os.path.basename(img_str)[:12])
   print(frac_lac_single_subject(nb.load(img_str.rstrip()).get_fdata()))
   print('Time taken for operation to run ', time.time() - t0 , 'seconds') 

elif choice == '2':
   print('Step 1: Enter paths to the folder containing the Images')

   print('Enter GBM Path: ')
   path_gbm = input() + '/'
   print('Enter LGG Path: ')
   path_lgg = input() + '/'

   path_list = [path_gbm, path_lgg]
   print('\n')
   print('*' * 50)
   
   inputfile_list = []
   for root_directory in path_list:      
      #List of GBM Subjects
      # root_directory = path_gbm
      file_pattern = '_GlistrBoost_ManuallyCorrected.nii.gz' #query

      file_paths = []

      # Traverse through the root directory and its subdirectories
      for root, dirs, files in os.walk(root_directory):
         for file in files:
            if file.endswith(file_pattern):
                  file_paths.append(os.path.join(root, file))

      if root_directory == path_gbm:
         output_file = 'gbm_input.txt'
      elif root_directory == path_lgg:
         output_file = 'lgg_input.txt'

      inputfile_list.append(output_file)

      with open(output_file, 'w') as f:
         for file_path in file_paths:
            f.write(file_path + '\n')

      print(f"{len(file_paths)} file paths were written to {output_file}")

   print('*' * 50, '\n')

   # inputfile_list = ['gbm_input.txt', 'lgg_input.txt']
   for sub_list in inputfile_list:
      print(f"Processing files in '{sub_list}'")      
      container = []
      with open(sub_list) as file:      #text file containing filenames
         for fl in file:
            try:
               container.append(frac_lac(nb.load(fl.rstrip()).get_fdata()))
            except:
               continue
      x = np.array(container)      
      case_id = []

      with open(sub_list) as file:
         for fl in file:
            case_id.append(os.path.basename(fl)[:12])    #Number is the position of subject code in that path filename
      print('Processed! Example of Case ID: ',case_id[0]) 

      if sub_list == 'gbm_input.txt':
         gbm_fract_dim = pd.DataFrame({'ID':case_id,
                        'ncr_net_meanlac': x[:, 0], 'ncr_net_medlac': x[:, 1],
                        'ed_meanlac': x[:, 2], 'ed_medlac':x[:, 3],
                        'et_meanlac': x[:, 4], 'et_medlac': x[:, 5]})
         gbm_fract_dim.to_csv('gbm_frac_dim.csv', index = False)

      elif sub_list == 'lgg_input.txt':
         lgg_fract_dim = pd.DataFrame({'ID':case_id,
                     'ncr_net_meanlac': x[:, 0], 'ncr_net_medlac': x[:, 1],
                     'ed_meanlac': x[:, 2], 'ed_medlac':x[:, 3],
                     'et_meanlac': x[:, 4], 'et_medlac': x[:, 5]})
         lgg_fract_dim.to_csv('lgg_frac_dim.csv', index = False)
   
   print('\n') 
   print('*' * 50)
   type_text('\n\nLacunarity calculation done! Writing results to Final_sheet_lac.csv!')

   df = pd.read_csv('gbm_frac_dim.csv')
   df2 = pd.read_csv('lgg_frac_dim.csv')

   merged = pd.concat([df, df2], axis = 0)
   merged.to_csv('Final_sheet_lac.csv', index = False)
   print('\n**** Final_sheet_lac.csv has been printed to the current directory! ****\n')
   
   files_to_remove = ['gbm_frac_dim.csv', 'lgg_frac_dim.csv']
   # print('Removing temporary files....')
   for file1 in files_to_remove:
    if os.path.exists(file1):         
        os.remove(file1)         
    else:
        print("Temporary file creation error! Suggested: Run code again or check the code for dependent packages.")
   # print('Temporary Files removed!')
   print('Operation completed in ', time.time() - t0 , ' seconds')
   type_text("Code Completed. Exiting.....")

elif choice not in ['1', '2']:
   print('Choose 1 or 2. Code has ended. Re-run code.')
   type_text("Code Completed. Exiting.....")
   

exit()