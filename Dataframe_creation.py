import subprocess
import os
import pandas as pd
import imageio
import imagej
import xlsxwriter
import sys
import openpyxl
import scyjava as sj
import csv
import numpy as np
from scipy.spatial import distance
import numpy as np
from skimage import io, filters
import matplotlib.pyplot as plt
import tifffile as tiff
from scipy.ndimage import binary_fill_holes
from skimage import io, filters, morphology, segmentation, color
from scipy.spatial.distance import cdist
from scipy.spatial import KDTree
from sklearn.neighbors import NearestNeighbors

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)


def process_tiff_file(tiff_file, output_tiff_file, counts, bins):
    images = io.imread(tiff_file)

    thresholds = [ostu_auto_thres(counts, bins) for slice in images]

    segmented_images = np.array(
        [(binary_fill_holes(slice > threshold)).astype(np.uint8) for slice, threshold in zip(images, thresholds)])
    segmented_images = (segmented_images * 255).astype(np.uint8)

    tiff.imwrite(output_tiff_file, segmented_images, dtype=np.uint8)

def test_write_permissions(directory):
    test_file_path = os.path.join(directory, "test_write_permissions.txt")
    try:
        with open(test_file_path, 'w') as test_file:
            test_file.write("This is a test file to check write permissions.")
        os.remove(test_file_path)
        print(f"Write permission to {directory} is confirmed.")
        return True
    except IOError:
        print(f"Write permission to {directory} is denied.")
        return False

def imagetocsv(file_path):
    # Convert the file path to use forward slashes
    file_path_corrected = file_path.replace("\\", "/")

    # Determine the output directory at the root of the project (outside .venv)
    image_name = os.path.splitext(os.path.basename(file_path_corrected))[0]
    output_directory = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"{image_name}-data")

    # Ensure the output directory exists
    os.makedirs(output_directory, exist_ok=True)

    # Check write permissions for the output directory
    if not test_write_permissions(output_directory):
        print(f"Cannot write to directory: {output_directory}")
        return

    # Initialize ImageJ in headless mode
    ij = imagej.init(['sc.fiji:fiji'], mode='headless')

    # Import necessary classes
    BF = sj.jimport('loci.plugins.BF')
    WindowManager = sj.jimport('ij.WindowManager')
    IJ = sj.jimport('ij.IJ')
    ResultsTable = sj.jimport('ij.measure.ResultsTable')

    try:
        # Open the image
        print(f"Opening image: {file_path_corrected}")
        imps = BF.openImagePlus(file_path_corrected)
        if not imps:
            print(f"Failed to open image: {file_path_corrected}")
            return

        for imp in imps:
            imp.show()

        # Run the macro to split channels
        split_macro = "run('Split Channels');"
        ij.py.run_macro(split_macro)

        # Get the list of all open windows
        window_titles = WindowManager.getImageTitles()
        # Identify the correct windows for each channel based on the observed titles
        channel_windows = {
            "micro": next((title for title in window_titles if "C3-" in title), None),
            "nu": next((title for title in window_titles if "C2-" in title), None),
            "cv": next((title for title in window_titles if "C1-" in title), None)
        }

        if None in channel_windows.values():
            print("Error: Could not find all channel windows.")
            return

        file_paths = {
            "micro": os.path.join(output_directory, f"{image_name}-micro.csv"),
            "nu": os.path.join(output_directory, f"{image_name}-nu.csv"),
            "cv": os.path.join(output_directory, f"{image_name}-cv.csv")
        }

        for filename, window_title in channel_windows.items():
            macro = f"""
            selectWindow("{window_title}");
            run("Z Project...", "projection=[Max Intensity]");
            getHistogram(values, counts, 256);
            for (i=0; i<values.length; i++) {{
                setResult("Value", i, values[i]);
                setResult("Count", i, counts[i]);
            }}
            updateResults();
            """
            ij.py.run_macro(macro)

            # Extract the histogram data from the ImageJ ResultsTable
            rt = ResultsTable.getResultsTable()
            values = rt.getColumn(0)
            counts = rt.getColumn(1)
            rt.reset()

            if values is None or counts is None:
                print(f"Error: No histogram data found for {filename}")
                continue

            # Combine bin edges and histogram counts into a list of tuples
            histogram_data = list(zip(values, counts))

            # Write histogram data to CSV file
            with open(file_paths[filename], 'w', newline='') as csvfile:
                csvwriter = csv.writer(csvfile)
                csvwriter.writerow(['', 'Bin', 'Count'])
                for bin_value, count in histogram_data:
                    csvwriter.writerow(['', f"{bin_value:.3f}", int(count)])
            print(f"Saved {filename}.csv to {file_paths[filename]}")

        print(f"Macro executed successfully. Files saved in {output_directory}")
        return file_paths
    except Exception as e:
        print(f"Error: {e}")
        return {}


def pdcsv_to_list(file_path):
    df = pd.read_csv(file_path)
    bin = df['Bin Start'].tolist()
    count = df['Count'].tolist()
    return bin, count

def csv_to_list(file_path):
    bin = []
    count = []
    with open(file_path, newline='') as csvfile:
        csvreader = csv.reader(csvfile)

        # Skip the header row if there is one
        next(csvreader)

        for row in csvreader:
            if len(row) >= 3:  # Ensure there are at least 3 columns
                bin.append(float(row[1]))  # 2nd column
                count.append(float(row[2]))  # 3rd column
    return bin, count

def man_thres(bin, count, abs_thres, correc_factor):
    total_auc = sum(count)
    auto_auc = 0.0
    man_auc = 0.0
    man_thres = 0.0
    counter = 1

    for i in range(len(bin) - 1, 0, -1):
        if bin[i] > abs_thres:
            counter += 1
        else:
            if abs(abs_thres - bin[i + 1]) > abs(abs_thres - bin[i]):
                counter += 1
            break
    for i in range(len(count) - 1, len(count) - counter, -1):
        auto_auc += count[i]
    auto_percent_auc = auto_auc / total_auc
    man_percent_auc = auto_percent_auc / correc_factor
    per_count = 1

    temp = 0.0
    per_count = 1
    for i in range(len(count) - 1, 0, -1):
        temp += count[i]
        if abs(man_percent_auc - (temp / total_auc)) > abs(man_percent_auc - ((temp + count[i - 1]) / total_auc)):
            per_count += 1
            man_auc = temp + count[i - 1]
        else:
            break

    man_thres = bin[(len(bin)) - per_count]

    return [total_auc, auto_auc, round(auto_percent_auc * 100, 5), round(man_percent_auc * 100, 5), man_auc,
            man_thres]

def data(mg, cv, nu, mgab, cvab, nuab):
    mgbin, mgcount = csv_to_list(mg)
    print(mgbin)
    print(mgcount)
    cvbin, cvcount = csv_to_list(cv)
    print(cvbin)
    print(cvcount)
    nubin, nucount = csv_to_list(nu)
    print(nubin)
    print(nucount)
    nudata = man_thres(nubin, nucount, nuab, 1.4)
    mgdata = man_thres(mgbin, mgcount, mgab, 1.227)
    print("Nucleus Manual Threshold for 1.4: " + str(nudata[5]))
    print("Mircogel Manual Threshold for 1.227: " + str(mgdata[5]))
    for i in range(8, 21):
        cf = i / 10.0
        mgdata = man_thres(mgbin, mgcount, mgab, cf)
        cudata = man_thres(cvbin, cvcount, cvab, cf)
        nudata = man_thres(nubin, nucount, nuab, cf)
        print("Correction Factor: " + str(cf))
        print("Mircogel Manual Threshold: " + str(mgdata[5]))
        print("Cell Volume Manual Threshold: " + str(cudata[5]))

#this funct takes only one channel images and turns them into csv files through the historgram while the other uses 3 channel data
def ch1imagetocsv(file_path):
    import os
    import csv
    import imagej
    import scyjava as sj

    # Convert the file path to use forward slashes
    file_path_corrected = file_path.replace("\\", "/")

    # Determine the output directory at the root of the project (outside .venv)
    image_name = os.path.splitext(os.path.basename(file_path_corrected))[0]
    output_directory = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"{image_name}-data")

    # Ensure the output directory exists
    os.makedirs(output_directory, exist_ok=True)

    # Check write permissions for the output directory
    def test_write_permissions(directory):
        try:
            testfile = os.path.join(directory, 'testfile')
            with open(testfile, 'w') as f:
                f.write('test')
            os.remove(testfile)
            return True
        except IOError:
            return False

    if not test_write_permissions(output_directory):
        print(f"Cannot write to directory: {output_directory}")
        return

    # Initialize ImageJ in headless mode
    ij = imagej.init(['sc.fiji:fiji'], mode='headless')

    # Import necessary classes
    BF = sj.jimport('loci.plugins.BF')
    WindowManager = sj.jimport('ij.WindowManager')
    IJ = sj.jimport('ij.IJ')
    ResultsTable = sj.jimport('ij.measure.ResultsTable')

    try:
        # Open the image
        print(f"Opening image: {file_path_corrected}")
        imps = BF.openImagePlus(file_path_corrected)
        if not imps:
            print(f"Failed to open image: {file_path_corrected}")
            return

        imp = imps[0]
        imp.show()

        file_path = os.path.join(output_directory, f"{image_name}-channel.csv")

        macro = f"""
        run("Z Project...", "projection=[Max Intensity]");
        getHistogram(values, counts, 256);
        for (i=0; i<values.length; i++) {{
            setResult("Value", i, values[i]);
            setResult("Count", i, counts[i]);
        }}
        updateResults();
        """
        ij.py.run_macro(macro)

        # Extract the histogram data from the ImageJ ResultsTable
        rt = ResultsTable.getResultsTable()
        values = rt.getColumn(0)
        counts = rt.getColumn(1)
        rt.reset()

        if values is None or counts is None:
            print(f"Error: No histogram data found for the image")
            return

        # Combine bin edges and histogram counts into a list of tuples
        histogram_data = list(zip(values, counts))

        # Write histogram data to CSV file
        with open(file_path, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(['', 'Bin', 'Count'])
            for bin_value, count in histogram_data:
                csvwriter.writerow(['', f"{bin_value:.3f}", int(count)])
        print(f"Saved channel.csv to {file_path}")

        print(f"Macro executed successfully. Files saved in {output_directory}")
        return file_path
    except Exception as e:
        print(f"Error: {e}")
        return {}

def ch1data(data, abs):
    bin, count = pdcsv_to_list(data)
    mgdata = man_thres(bin, count, abs, 1.5)
    print(data)
    for i in range(10, 21):
        cf = i / 10
        mgdata = man_thres(bin, count, abs, cf)
        print("Correction Factor: " + str(cf))
        print("Manual Threshold: " + str(mgdata[5]))

def process_csv_files(folder_path):
    # Iterate over all files and subfolders in the folder
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            # Check if the file is a CSV file
            if file.endswith(".csv"):
                file_path = os.path.join(root, file)

                try:
                    # Read the CSV file into a DataFrame, skipping the first three rows
                    df = pd.read_csv(file_path, skiprows=3)

                    # Save the modified DataFrame back to the CSV file
                    df.to_csv(file_path, index=False)
                    print(f"Processed file: {file_path}")

                except Exception as e:
                    print(f"Error processing file {file_path}: {e}")

def ostu_auto_thres(counts, bins):
    """
    Compute the Otsu's threshold using counts and bins of the histogram.
    """
    # Total number of pixels
    total = np.sum(counts)

    # Initialize variables
    sumB = 0.0
    wB = 0.0
    maximum = 0.0
    sum1 = np.dot(bins, counts)
    threshold = 0
    wF = 0.0

    for i in range(len(counts)):
        wB += counts[i]  # Weight Background
        if wB == 0:
            continue

        wF = total - wB  # Weight Foreground
        if wF == 0:
            break

        sumB += bins[i] * counts[i]
        mB = sumB / wB  # Mean Background
        mF = (sum1 - sumB) / wF  # Mean Foreground

        # Calculate Between Class Variance
        varBetween = wB * wF * (mB - mF) ** 2

        # Check if new maximum found
        if varBetween > maximum:
            threshold = bins[i]
            maximum = varBetween

    return threshold

I214 = r"C:\Research Files\UIC - Shin Lab\Single Cell Encaps\IMARIS_Images_sorting\Stiff\Image_2\Image_2_Cell_volume_1.4_Position_Detailed.csv"
I214V = r"C:\Research Files\UIC - Shin Lab\Single Cell Encaps\IMARIS_Images_sorting\Stiff\Image_2\Image_2_Cell_volume_1.4_Detailed.csv"
I217V = r"C:\Research Files\UIC - Shin Lab\Single Cell Encaps\IMARIS_Images_sorting\Stiff\Image_2\CV-manual_1.7_image 2_Detailed.csv"
I217 = r"C:\Research Files\UIC - Shin Lab\Single Cell Encaps\IMARIS_Images_sorting\Stiff\Image_2\Image_2_Cell_1.7_manual_position_Detailed.csv"
I215 = r"C:\Research Files\UIC - Shin Lab\Single Cell Encaps\IMARIS_Images_sorting\Stiff\Image_2\Image_2_Cell_1.5_manual_position_Detailed.csv"
I216 = r"C:\Research Files\UIC - Shin Lab\Single Cell Encaps\IMARIS_Images_sorting\Stiff\Image_2\Image_2_Cell_1.6_Position_Detailed.csv"
I218 =r"C:\Research Files\UIC - Shin Lab\Single Cell Encaps\IMARIS_Images_sorting\Stiff\Image_2\Image_2_Cell_1.8_manual_position_Detailed.csv"
I220 = r"C:\Research Files\UIC - Shin Lab\Single Cell Encaps\IMARIS_Images_sorting\Stiff\Image_2\Image_2_Cell_2.0_manual_position_Detailed.csv"

I215V = r"C:\Research Files\UIC - Shin Lab\Single Cell Encaps\IMARIS_Images_sorting\Stiff\Image_2\Cell_Volume-Manual_1.5_correction_Detailed.csv"
I216V = r"C:\Research Files\UIC - Shin Lab\Single Cell Encaps\IMARIS_Images_sorting\Stiff\Image_2\Image_2_Cell_volume_1.6-manual-volume_Detailed.csv"
I218V = r"C:\Research Files\UIC - Shin Lab\Single Cell Encaps\IMARIS_Images_sorting\Stiff\Image_2\CV-manual_1.8_image_2_Detailed.csv"
I220V = r"C:\Research Files\UIC - Shin Lab\Single Cell Encaps\IMARIS_Images_sorting\Stiff\Image_2\Cell_volume-Manual_2.0_Correction_Detailed.csv"



def read_csv(file1, file2):
    df = pd.read_csv(file1)
    vdf = pd.read_csv(file2)

    vdf.drop(columns=['Unit', 'Category', 'Time',], inplace=True)
    df.drop(columns=['Unit', 'Category', 'Time', 'Collection'], inplace=True)
    df = df.loc[:, ['ID', 'Position X', 'Position Y', 'Position Z']]
    df['ID'] = df['ID'] - df['ID'][0] + 1
    vdf = vdf.loc[:, ['ID', 'Volume']]

    vdf['ID'] = vdf['ID'] - vdf['ID'][0] + 1
    new_data = pd.merge(vdf, df, on='ID')

    return new_data
I2CV14 = read_csv(I214, I214V)
I2CV15 = read_csv(I215, I215V)
I2CV16 = read_csv(I216, I216V)
I2CV17 = read_csv(I217, I217V)
I2CV18 = read_csv(I218, I218V)
I2CV20 = read_csv(I220, I220V)
max_l = max(I2CV14.shape, I2CV15.shape, I2CV16.shape, I2CV17.shape, I2CV17.shape, I2CV18.shape, I2CV20.shape)[0]
print(I2CV14.shape, I2CV15.shape, I2CV16.shape, I2CV17.shape, I2CV17.shape, I2CV18.shape, I2CV20.shape)
def padding(df, max_l=max_l):
    temp = np.zeros([max_l - df.shape[0], df.shape[1]])
   # temp = temp * np.NaN
#   for i in range(1, (max_l - df.shape[0]) + 1):
#        temp[i - 1, 0] = df.shape[0] + i

    temp = {
        'ID': temp[:, 0],
        'Volume': temp[:, 1],
        'Position X': temp[:, 2],
        'Position Y': temp[:, 3],
        'Position Z': temp[:, 4],
    }
    temp = pd.DataFrame(temp)
    temp =  pd.concat([df, temp], axis=0).reset_index(drop=True)
    return temp

I2CV14 = padding(I2CV14)
I2CV15 = padding(I2CV15)
print(I2CV14)
I2CV16 = padding(I2CV16)
I2CV17 = padding(I2CV17)
I2CV18 = padding(I2CV18)
I2CV20 = padding(I2CV20)

tree = KDTree(I2CV14[['Position X', 'Position Y', 'Position Z']])
tree1 = KDTree(I2CV15[['Position X', 'Position Y', 'Position Z']])
tree2 = KDTree(I2CV16[['Position X', 'Position Y', 'Position Z']])
tree3 = KDTree(I2CV17[['Position X', 'Position Y', 'Position Z']])
tree4 = KDTree(I2CV18[['Position X', 'Position Y', 'Position Z']])

dis1, ind1 = tree.query(I2CV15[['Position X', 'Position Y', 'Position Z']], k=1)
dis2, ind2 = tree.query(I2CV16[['Position X', 'Position Y', 'Position Z']], k=1)
dis3, ind3 = tree.query(I2CV17[['Position X', 'Position Y', 'Position Z']], k=1)
dis4, ind4 = tree.query(I2CV18[['Position X', 'Position Y', 'Position Z']], k=1)
dis5, ind5 = tree.query(I2CV20[['Position X', 'Position Y', 'Position Z']], k=1)

m_1415 = pd.DataFrame({
   'CV14': I2CV14.iloc[ind1]['ID'],
   'CV15': I2CV15['ID'].values,
}).sort_values(by='CV14').reset_index(drop=True)

#print(I2CV14.sort_values(by='Position X')) # 16  17.0   3417.640    205.7850   645.48500    58.47470
print(I2CV15.sort_values(by='Position X')) # 16  17.0   3417.640    205.7850   645.48500    58.47470

print(I2CV14.loc[I2CV14['ID'] == 17])
print(I2CV15.loc[I2CV15['ID'] == 14]) # 1817.34     206.104     640.825     64.1261
print(I2CV15.loc[I2CV15['ID'] == 17])
print(m_1415.sort_values(by='CV14').reset_index(drop=True))
















I2CV14 = I2CV14.sort_values(by='Position X', ascending=False).reset_index(drop=True)
I2CV15 = I2CV15.sort_values(by='Position X', ascending=False).reset_index(drop=True)
I2CV16 = I2CV16.sort_values(by='Position X', ascending=False).reset_index(drop=True)
I2CV17 = I2CV17.sort_values(by='Position X', ascending=False).reset_index(drop=True)
I2CV18 = I2CV18.sort_values(by='Position X', ascending=False).reset_index(drop=True)
I2CV20 = I2CV20.sort_values(by='Position X', ascending=False).reset_index(drop=True)

CV14 = pd.concat([I2CV14['ID'], I2CV14['Position X']], axis=1)
CV15 = pd.concat([I2CV15['ID'], I2CV15['Position X']], axis=1)
CV16 = pd.concat([I2CV16['ID'], I2CV16['Position X']], axis=1)
CV17 = pd.concat([I2CV17['ID'], I2CV17['Position X']], axis=1)
CV18 = pd.concat([I2CV18['ID'], I2CV18['Position X']], axis=1)
CV20 = pd.concat([I2CV20['ID'], I2CV20['Position X']], axis=1)

