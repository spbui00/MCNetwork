import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def read_data(file_path):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            values = line.strip().split()
            data.append([float(v) for v in values])
    return np.array(data)

def analyze_data(data):
    current_values = data[:, -2]
    uncertainty_values = data[:, -1]

    stats = {
        'mean_current': np.mean(current_values),
        'median_current': np.median(current_values),
        'std_current': np.std(current_values),
        'var_current': np.var(current_values),
        'min_current': np.min(current_values),
        'max_current': np.max(current_values),
        'mean_uncertainty': np.mean(uncertainty_values),
        'median_uncertainty': np.median(uncertainty_values),
        'std_uncertainty': np.std(uncertainty_values),
        'var_uncertainty': np.var(uncertainty_values),
        'min_uncertainty': np.min(uncertainty_values),
        'max_uncertainty': np.max(uncertainty_values),
    }

    print("Current Statistics:")
    print(f"Mean: {stats['mean_current']:.2f} nA")
    print(f"Median: {stats['median_current']:.2f} nA")
    print(f"Standard Deviation: {stats['std_current']:.2f} nA")
    print(f"Variance: {stats['var_current']:.2f} (nA)^2")
    print(f"Min: {stats['min_current']:.2f} nA")
    print(f"Max: {stats['max_current']:.2f} nA")
    print()
    print("Uncertainty Statistics:")
    print(f"Mean: {stats['mean_uncertainty']:.2f} nA")
    print(f"Median: {stats['median_uncertainty']:.2f} nA")
    print(f"Standard Deviation: {stats['std_uncertainty']:.2f} nA")
    print(f"Variance: {stats['var_uncertainty']:.2f} (nA)^2")
    print(f"Min: {stats['min_uncertainty']:.2f} nA")
    print(f"Max: {stats['max_uncertainty']:.2f} nA")

    return current_values, uncertainty_values, stats

def plot_histogram(data, title, xlabel, ylabel, bins=30):
    plt.figure()
    plt.hist(data, bins=bins, edgecolor='black')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.show()

def plot_error_bars(current_values, uncertainty_values):
    plt.figure()
    plt.errorbar(range(len(current_values)), current_values, yerr=uncertainty_values, fmt='o', ecolor='r', capsize=5, label='Current with Uncertainty')
    plt.title('Output Currents with Uncertainties')
    plt.xlabel('Sample Index')
    plt.ylabel('Current (nA)')
    plt.legend()
    plt.grid(True)
    plt.show()

def main():
    uncertainties_file = '../data/samples/1e5/IO_uncertainties.dat'

    data = read_data(uncertainties_file)

    current_values, uncertainty_values, stats = analyze_data(data)

    plot_histogram(uncertainty_values, 'Histogram of Uncertainties', 'Uncertainty (nA)', 'Frequency')

if __name__ == "__main__":
    main()