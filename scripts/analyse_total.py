import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def analyze_io_file(file_path):
    # Read the data from the file
    data = pd.read_csv(file_path, comment='#', header=None, delim_whitespace=True)

    # Extract inputs and outputs
    inputs = data.iloc[:, :-1]
    outputs = data.iloc[:, -1]

    # Calculate statistical measures
    mean = outputs.mean()
    std_dev = outputs.std()
    variance = outputs.var()
    data_range = outputs.max() - outputs.min()
    median = outputs.median()
    iqr = outputs.quantile(0.75) - outputs.quantile(0.25)

    # Print the results
    print(f"Mean: {mean}")
    print(f"Standard Deviation: {std_dev}")
    print(f"Variance: {variance}")
    print(f"Range: {data_range}")
    print(f"Median: {median}")
    print(f"Interquartile Range (IQR): {iqr}")

    # Plotting
    plt.figure(figsize=(14, 10))

    # Histogram
    plt.subplot(2, 2, 1)
    plt.hist(outputs, bins=30, edgecolor='black')
    plt.title('Histogram of Outputs')
    plt.xlabel('Output Value')
    plt.ylabel('Frequency')

    # Box Plot
    plt.subplot(2, 2, 2)
    plt.boxplot(outputs, vert=False)
    plt.title('Box Plot of Outputs')
    plt.xlabel('Output Value')

    # Hexbin Plot of Inputs vs Output
    plt.subplot(2, 2, 3)
    for i in range(inputs.shape[1]):
        plt.hexbin(inputs.iloc[:, i], outputs, gridsize=50, cmap='viridis', alpha=0.5)
    plt.title('Hexbin Plot of Inputs vs Output')
    plt.xlabel('Input Value')
    plt.ylabel('Output Value')

    # Show the first set of plots
    plt.tight_layout()
    plt.show()

    # Pair Plot (using a subset of the data for better visualization)
    subset_data = pd.concat([inputs, outputs], axis=1).sample(n=1000, random_state=1)
    sns.pairplot(subset_data)
    plt.suptitle('Pair Plot of Inputs and Output', y=1.02)
    plt.show()

def analyze_with_same_inputs(main_file_path, same_input_file_path):
    # Read the data from the files
    random_data = pd.read_csv(main_file_path, comment='#', header=None, delim_whitespace=True)
    same_input_data = pd.read_csv(same_input_file_path, comment='#', header=None, delim_whitespace=True)

    # Extract outputs
    random_outputs = random_data.iloc[:, -1]
    same_input_outputs = same_input_data.iloc[:, -1]

    # Calculate statistical measures for random outputs
    random_mean = random_outputs.mean()
    random_std_dev = random_outputs.std()
    random_variance = random_outputs.var()
    random_range = random_outputs.max() - random_outputs.min()
    random_median = random_outputs.median()
    random_iqr = random_outputs.quantile(0.75) - random_outputs.quantile(0.25)

    # Calculate statistical measures for same input outputs
    same_input_mean = same_input_outputs.mean()
    same_input_std_dev = same_input_outputs.std()
    same_input_variance = same_input_outputs.var()
    same_input_range = same_input_outputs.max() - same_input_outputs.min()
    same_input_median = same_input_outputs.median()
    same_input_iqr = same_input_outputs.quantile(0.75) - same_input_outputs.quantile(0.25)

    # Print the results
    print("Random Outputs Statistical Measures:")
    print(f"Mean: {random_mean}")
    print(f"Standard Deviation: {random_std_dev}")
    print(f"Variance: {random_variance}")
    print(f"Range: {random_range}")
    print(f"Median: {random_median}")
    print(f"Interquartile Range (IQR): {random_iqr}")
    print()
    print("Same Input Outputs Statistical Measures:")
    print(f"Mean: {same_input_mean}")
    print(f"Standard Deviation: {same_input_std_dev}")
    print(f"Variance: {same_input_variance}")
    print(f"Range: {same_input_range}")
    print(f"Median: {same_input_median}")
    print(f"Interquartile Range (IQR): {same_input_iqr}")

    # Histogram for Same Input Outputs
    plt.subplot(2, 2, 2)
    plt.hist(same_input_outputs, bins=30, edgecolor='black')
    plt.title('Histogram of Same Input Outputs')
    plt.xlabel('Output Value')
    plt.ylabel('Frequency')

    # Box Plot for Random Outputs
    plt.subplot(2, 2, 3)
    plt.boxplot(random_outputs, vert=False)
    plt.title('Box Plot of Random Outputs')
    plt.xlabel('Output Value')

    # Box Plot for Same Input Outputs
    plt.subplot(2, 2, 4)
    plt.boxplot(same_input_outputs, vert=False)
    plt.title('Box Plot of Same Input Outputs')
    plt.xlabel('Output Value')

    # Show the plots
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze output files.")
    parser.add_argument("--total", type=str, help="Path to the output file to analyze.")
    parser.add_argument("--consistency", action='store_true', help="Flag to analyze consistency.")
    parser.add_argument("--mainFile", type=str, help="Path to the main file with all data.")
    # parser.add_argument("--sameInputFile", type=str, help="Path to the file with same inputs.")

    args = parser.parse_args()

    if args.total:
        analyze_io_file(args.total)
    elif args.consistency and args.mainFile and args.sameInputFile:
        analyze_with_same_inputs(args.mainFile, args.sameInputFile)
    else:
        print("Invalid arguments. Use --total to analyze a single output file or --consistency with --mainFile and --sameInputFile to analyze consistency.")