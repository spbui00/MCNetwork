import numpy as np
import matplotlib.pyplot as plt
import argparse

def analyze_io_file(file_path):
    # Read the data from the file
    data = np.loadtxt(file_path, comments='#')

    # Extract the outputs (last column)
    outputs = data[:, -1]

    # Calculate statistical measures
    mean = np.mean(outputs)
    std_dev = np.std(outputs)
    variance = np.var(outputs)
    cv = std_dev / mean if mean != 0 else float('inf')
    data_range = np.ptp(outputs)
    median = np.median(outputs)

    # Print the results
    print(f"Mean: {mean}")
    print(f"Standard Deviation: {std_dev}")
    print(f"Variance: {variance}")
    print(f"Coefficient of Variation: {cv}")
    print(f"Range: {data_range}")
    print(f"Median: {median}")

    # Plotting
    plt.figure(figsize=(12, 8))

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

    # Time Series Plot
    plt.subplot(2, 2, 3)
    plt.plot(outputs, marker='o', linestyle='-', markersize=2)
    plt.title('Time Series of Outputs')
    plt.xlabel('Sample Index')
    plt.ylabel('Output Value')

    # Show plots
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Analyze IO.dat file.')
    parser.add_argument('file_path', type=str, help='Path to the IO.dat file')
    args = parser.parse_args()

    analyze_io_file(args.file_path)