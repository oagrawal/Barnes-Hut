import subprocess
import matplotlib.pyplot as plt
import numpy as np
import time
import re
import os

def run_command(command):
    # Run the command and extract the execution time
    print("  Running command...")
    process = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    output = process.stdout.strip()
    
    try:
        # Directly parse the output as a float
        time_value = float(output)
        print(f"  Execution time: {time_value:.6f} seconds")
        return time_value
    except ValueError:
        print(f"  Failed to extract time from output: '{output}'")
        print(f"  STDERR: '{process.stderr}'")
        return None

def warmup():
    # Run a warm-up command using nb-10.txt
    print("Warming up processors using input/nb-10.txt...")
    warmup_command = "HWLOC_COMPONENTS=-gl mpirun -n 2 --oversubscribe ./nbody -i input/nb-10.txt -o output/output.txt -s 20 -t 0.35 -d 0.005"
    subprocess.run(warmup_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print("Warmup complete")

def ensure_output_directory():
    # Make sure the output directory exists
    os.makedirs("output", exist_ok=True)

def plot_speedup(results, worker_counts):
    if results["sequential"] is None:
        return
    
    workers = []
    speedups = []
    
    for w in worker_counts:
        if w in results["parallel"] and results["parallel"][w] is not None:
            workers.append(w)
            speedups.append(results["sequential"] / results["parallel"][w])
    
    if not workers:
        print("No valid parallel results, skipping plot")
        return
    
    plt.figure(figsize=(10, 6))
    plt.plot(workers, speedups, marker='o', linestyle='-', linewidth=2, markersize=10, color='blue')
    plt.title("Speedup vs. Number of Workers for input/nb-100000.txt", fontsize=14)
    plt.xlabel("Number of Workers", fontsize=12)
    plt.ylabel("Speedup (Sequential Time / Parallel Time)", fontsize=12)
    plt.grid(True)
    
    # Add a horizontal line at speedup = 1 (no speedup)
    # plt.axhline(y=1, color='r', linestyle='--', label='No Speedup')
    
    # Add linear speedup reference line
    max_w = max(workers)
    # plt.plot([1, max_w], [1, max_w], 'g--', label='Ideal Linear Speedup')
    
    plt.legend(fontsize=10)
    plt.xticks(workers)
    
    # Add data point annotations
    for i, (w, s) in enumerate(zip(workers, speedups)):
        plt.annotate(f"{s:.2f}x", 
                    (w, s), 
                    textcoords="offset points",
                    xytext=(0, 10), 
                    ha='center',
                    fontsize=10)
    
    plt.tight_layout()
    
    # Save the plot
    plot_filename = "speedup_nb-100000.png"
    plt.savefig(plot_filename, dpi=300)
    print(f"Saved speedup plot to {plot_filename}")

def main():
    # Start timing the entire script
    start_time = time.time()
    
    # Ensure output directory exists
    ensure_output_directory()
    
    test_case = "input/nb-100000.txt"
    worker_counts = [8]  # Just one worker count for now
    
    results = {
        "sequential": None,
        "parallel": {}
    }
    
    print(f"\nProcessing test case: {test_case}")
    
    # Warm up
    warmup()
    
    # Run sequential implementation
    print("\nRunning sequential implementation...")
    seq_command = f"HWLOC_COMPONENTS=-gl mpirun -n 1 --oversubscribe ./nbody -i {test_case} -o output/output.txt -s 20 -t 0.35 -d 0.005 -S"
    results["sequential"] = run_command(seq_command)
    
    if results["sequential"] is None:
        print(f"Failed to run sequential implementation for {test_case}")
        return
    
    print(f"Sequential time: {results['sequential']:.6f} seconds")
    
    # Run parallel implementation with specified worker counts
    for workers in worker_counts:
        # Warm up before each test
        warmup()
        
        print(f"\nRunning parallel implementation with {workers} workers...")
        par_command = f"HWLOC_COMPONENTS=-gl mpirun -n {workers} --oversubscribe ./nbody -i {test_case} -o output/output.txt -s 20 -t 0.35 -d 0.005"
        results["parallel"][workers] = run_command(par_command)
        
        if results["parallel"][workers] is None:
            print(f"Failed to run parallel implementation with {workers} workers")
            continue
        
        speedup = results["sequential"] / results["parallel"][workers]
        print(f"Parallel time with {workers} workers: {results['parallel'][workers]:.6f} seconds, Speedup: {speedup:.2f}x")
    
    # Plot speedup graph
    plot_speedup(results, worker_counts)
    
    # Display detailed results in terminal
    if results["sequential"] is not None:
        print(f"\nDetailed Results for {test_case}:")
        print(f"{'Workers':<10} {'Time (s)':<15} {'Speedup':<10}")
        print("-" * 35)
        print(f"{'Sequential':<10} {results['sequential']:<15.6f} {1.0:<10.2f}")
        
        for w in worker_counts:
            if w in results["parallel"] and results["parallel"][w] is not None:
                speedup = results["sequential"] / results["parallel"][w]
                print(f"{w:<10} {results['parallel'][w]:<15.6f} {speedup:<10.2f}")
    
    # Save the results to a file
    with open("performance_results_nb100000.txt", "w") as f:
        f.write(f"Test Case: {test_case}\n")
        if results["sequential"] is not None:
            f.write(f"Sequential Time: {results['sequential']:.6f} seconds\n")
            f.write("Parallel Results:\n")
            f.write(f"{'Workers':<10} {'Time (s)':<15} {'Speedup':<10}\n")
            f.write("-" * 35 + "\n")
            
            for w in worker_counts:
                if w in results["parallel"] and results["parallel"][w] is not None:
                    speedup = results["sequential"] / results["parallel"][w]
                    f.write(f"{w:<10} {results['parallel'][w]:<15.6f} {speedup:<10.2f}\n")
        else:
            f.write("Sequential implementation failed, no results available.\n")
    
    # Calculate and display total script execution time
    total_time = time.time() - start_time
    hours, remainder = divmod(total_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    
    print("\n" + "="*50)
    print(f"TOTAL SCRIPT EXECUTION TIME: {int(hours):02d}:{int(minutes):02d}:{seconds:.2f}")
    print(f"                             ({total_time:.2f} seconds)")
    print("="*50)
    
    # Also add the timing information to the results file
    with open("performance_results_nb100000.txt", "a") as f:
        f.write(f"\nTotal script execution time: {int(hours):02d}:{int(minutes):02d}:{seconds:.2f} ({total_time:.2f} seconds)\n")
    
    print("\nAll results saved to performance_results_nb100000.txt")

if __name__ == "__main__":
    main()
