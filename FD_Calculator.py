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


#Fractal Dimension Calculation function
def fract(x):                                             
   x = x.astype(np.uint8)
   y = x      
   bin_ = np.zeros((x.shape[0], x.shape[1]))    
   bin_[x!=0] = 1
   bin_ = bin_.astype(np.uint8)    
   contours = cv2.findContours(image = bin_, mode = cv2.RETR_TREE, method = cv2.CHAIN_APPROX_SIMPLE)[0]
   contours = sorted(contours, key = cv2.contourArea, reverse = True)

   if len(contours) == 0:
      return 0

   c_0 = contours[0]                                   # Get the 4 points of the bounding rectangle
   xx, yy, w, h = cv2.boundingRect(c_0)    
   X  =  y[yy:yy+h, xx:xx+w]    
      
   XX = np.zeros((X.shape[0], X.shape[1]))       
   XX[X==1] = 1    
   ncrnet =  fract_dim(XX)    

   XX = np.zeros((X.shape[0], X.shape[1]))    
   XX[X==2] = 1    
   ed =  fract_dim(XX)    

   XX = np.zeros((X.shape[0], X.shape[1]))    
   XX[X==4] = 1    
   et =  fract_dim(XX)    

   return np.array([ncrnet, ed, et])

def final_fd_single_search(X):
    
   X = X.astype(np.uint8)
   #print('pre removing empty slice',X.shape)

   X = np.delete(X, [i for i in range(0, 155) if sum(X[:,:,i].ravel()) < 10], 2)
   #print('post removing empty slice',X.shape)
   
   X = [X[:,:,i] for i in range(0, X.shape[-1])]

   pool = mp.Pool(processes = 3) 
   fd = np.array(pool.map(fract, X))
   
   n_mean =  np.mean(fd[:,0])
   n_med =  np.median(fd[:,0]) 
   
   ed_mean =  np.mean(fd[:,1])
   ed_med =  np.median(fd[:,1])   
   
   et_mean =  np.mean(fd[:,2])
   et_med =  np.median(fd[:,2])

   result = np.array([n_mean, n_med, ed_mean, ed_med, et_mean, et_med])
   print('Results: ')
   print('Non enhancing FD MEAN: ', result[0], '\nEdema FD MEAN: ', result[2], '\nEnhancing FD MEAN: ', result[4])
   return type_text("Code Completed")

def final_fd_(X):
    
   X = X.astype(np.uint8)
   #print('pre removing empty slice',X.shape)

   X = np.delete(X, [i for i in range(0, 155) if sum(X[:,:,i].ravel()) < 10], 2)
   #print('post removing empty slice',X.shape)
   
   X = [X[:,:,i] for i in range(0, X.shape[-1])]

   pool = mp.Pool(processes = 3) 
   fd = np.array(pool.map(fract, X))
   
   n_mean =  np.mean(fd[:,0])
   n_med =  np.median(fd[:,0]) 
   
   ed_mean =  np.mean(fd[:,1])
   ed_med =  np.median(fd[:,1])   
   
   et_mean =  np.mean(fd[:,2])
   et_med =  np.median(fd[:,2])

   return  np.array([n_mean, n_med, ed_mean, ed_med, et_mean, et_med])

type_text('Done!')
print('')
print('\n*******************')
print('')

print('1. Check code integrity for a single subject.  \n2. Calculate batch fractal dimension. \n\nType 1 or 2 to choose.')

choice = str(input())

if choice == '1':
   print('\n')
   print('*' * 50)
   print('File path for the image file (Filename ends with: _Glistrboost_Manually_Corrected.nii.gz): ')
   # print('Example: /home/username/directory/file_glistrboost_manually_corrected.nii.gz (or ".nii" file)')
   img_str = str(input())
   print(os.path.basename(img_str)[:12] , final_fd_single_search(nb.load(img_str.rstrip()).get_fdata()))
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
               container.append(final_fd_(nb.load(fl.rstrip()).get_fdata()))
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
                        'ncr_net_meanfd': x[:, 0], 'ncr_net_medfd': x[:, 1],
                        'ed_meanfd': x[:, 2], 'ed_medfd':x[:, 3],
                        'et_meanfd': x[:, 4], 'et_medfd': x[:, 5]})
         gbm_fract_dim.to_csv('gbm_frac_dim.csv', index = False)

      elif sub_list == 'lgg_input.txt':
         lgg_fract_dim = pd.DataFrame({'ID':case_id,
                     'ncr_net_meanfd': x[:, 0], 'ncr_net_medfd': x[:, 1],
                     'ed_meanfd': x[:, 2], 'ed_medfd':x[:, 3],
                     'et_meanfd': x[:, 4], 'et_medfd': x[:, 5]})
         lgg_fract_dim.to_csv('lgg_frac_dim.csv', index = False)
   
   print('\n') 
   print('*' * 50)
   type_text('\n\nFD calculation done! Writing results to Final_sheet.csv!')

   df = pd.read_csv('gbm_frac_dim.csv')
   df2 = pd.read_csv('lgg_frac_dim.csv')

   merged = pd.concat([df, df2], axis = 0)
   merged.to_csv('Final_sheet.csv', index = False)
   print('\n**** Final_sheet.csv has been printed to the current directory! ****\n')
   
   files_to_remove = ['gbm_frac_dim.csv', 'lgg_frac_dim.csv']
   # print('Removing temporary files....')
   for file1 in files_to_remove:
    if os.path.exists(file1):         
        os.remove(file1)         
    else:
        print("Temporary file creation error! Suggested: Run code again or check the code for dependent packages.")
   # print('Temporary Files removed!')
   print('Operation completed in ', time.time() - t0 , ' seconds') 

elif choice not in ['1', '2']:
   print('Choose 1 or 2. Code has ended. Re-run code.')
    
exit()