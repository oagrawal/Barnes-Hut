import subprocess
import matplotlib.pyplot as plt
import numpy as np
import time
import re
import os

def run_command(command):
    # Run the command and extract the execution time
    times = []
    for i in range(3):  # Run 3 times
        print(f"  Run {i+1}/3...")
        process = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        output = process.stdout.strip()
        
        try:
            # Directly parse the output as a float
            time_value = float(output)
            times.append(time_value)
            print(f"  Execution time: {time_value:.6f} seconds")
        except ValueError:
            print(f"  Failed to extract time from output: '{output}'")
            print(f"  STDERR: '{process.stderr}'")
            times.append(None)
        
        # Small delay between runs
        time.sleep(1)
    
    # Calculate the average time
    valid_times = [t for t in times if t is not None]
    if valid_times:
        avg_time = sum(valid_times) / len(valid_times)
        print(f"  Average time: {avg_time:.6f} seconds")
        return avg_time
    else:
        print("  No valid times recorded")
        return None

def warmup(test_case):
    # Run a warm-up command to prepare processors
    print(f"Warming up processors using {test_case}...")
    warmup_command = f"HWLOC_COMPONENTS=-gl mpirun -n 2 --oversubscribe ./nbody -i {test_case} -o output/output.txt -s 100 -t 0.35 -d 0.005"
    subprocess.run(warmup_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print("Warmup complete")

def ensure_output_directory():
    # Make sure the output directory exists
    os.makedirs("output", exist_ok=True)

def plot_speedup(results, test_case, worker_counts):
    if results[test_case]["sequential"] is None:
        return
    
    workers = []
    speedups = []
    
    for w in worker_counts:
        if w in results[test_case]["parallel"] and results[test_case]["parallel"][w] is not None:
            workers.append(w)
            speedups.append(results[test_case]["sequential"] / results[test_case]["parallel"][w])
    
    if not workers:
        print(f"No valid parallel results for {test_case}, skipping plot")
        return
    
    plt.figure(figsize=(12, 7))
    plt.plot(workers, speedups, marker='o', linestyle='-', linewidth=2, markersize=8)
    plt.title(f"Speedup vs. Number of Workers for {test_case}", fontsize=14)
    plt.xlabel("Number of Workers", fontsize=12)
    plt.ylabel("Speedup (Sequential Time / Parallel Time)", fontsize=12)
    plt.grid(True)
    
    # Add a horizontal line at speedup = 1 (no speedup)
    # plt.axhline(y=1, color='r', linestyle='--', label='No Speedup')
    
    # Add linear speedup reference line
    max_w = max(workers)
    # plt.plot([1, max_w], [1, max_w], 'g--', label='Ideal Linear Speedup')
    
    plt.legend(fontsize=10)
    plt.xticks(workers, rotation=45 if max_w > 20 else 0)
    
    # Add data point annotations
    for i, (w, s) in enumerate(zip(workers, speedups)):
        plt.annotate(f"{s:.2f}x", 
                    (w, s), 
                    textcoords="offset points",
                    xytext=(0, 10), 
                    ha='center',
                    fontsize=9)
    
    plt.tight_layout()
    
    # Save the plot
    plot_filename = f"speedup_{test_case.split('/')[-1].replace('.txt', '')}.png"
    plt.savefig(plot_filename, dpi=300)
    print(f"Saved speedup plot to {plot_filename}")

def main():
    # Ensure output directory exists
    ensure_output_directory()
    
    test_cases = ["input/nb-10.txt", "input/nb-100.txt"]
    
    # Define specific worker counts for each test case
    worker_counts = {
        "input/nb-10.txt": [1, 2, 4, 6, 8, 10],
        "input/nb-100.txt": [1, 2, 4, 6, 8, 10, 12, 16, 20, 25, 20, 25, 30, 35, 40, 45, 50]
    }
    
    results = {}
    
    for test_case in test_cases:
        print(f"\nProcessing test case: {test_case}")
        results[test_case] = {
            "sequential": 0,
            "parallel": {}
        }
        
        # Warm up
        warmup(test_case)
        
        # Run sequential implementation
        print("\nRunning sequential implementation...")
        seq_command = f"HWLOC_COMPONENTS=-gl mpirun -n 1 --oversubscribe ./nbody -i {test_case} -o output/output.txt -s 1000 -t 0.35 -d 0.005 -S"
        results[test_case]["sequential"] = run_command(seq_command)
        
        if results[test_case]["sequential"] is None:
            print(f"Failed to run sequential implementation for {test_case}")
            continue
        
        print(f"Sequential time for {test_case}: {results[test_case]['sequential']:.6f} seconds")
        
        # Run parallel implementation with specified worker counts
        for workers in worker_counts[test_case]:
            # Skip if we've already tested this worker count (for duplicates in nb-100 list)
            if workers in results[test_case]["parallel"]:
                print(f"\nSkipping already tested worker count: {workers}")
                continue
                
            # Warm up before each test
            warmup(test_case)
            
            print(f"\nRunning parallel implementation with {workers} workers...")
            par_command = f"HWLOC_COMPONENTS=-gl mpirun -n {workers} --oversubscribe ./nbody -i {test_case} -o output/output.txt -s 1000 -t 0.35 -d 0.005"
            results[test_case]["parallel"][workers] = run_command(par_command)
            
            if results[test_case]["parallel"][workers] is None:
                print(f"Failed to run parallel implementation with {workers} workers for {test_case}")
                continue
            
            speedup = results[test_case]["sequential"] / results[test_case]["parallel"][workers]
            print(f"Parallel time with {workers} workers: {results[test_case]['parallel'][workers]:.6f} seconds, Speedup: {speedup:.2f}x")
    
    # Plot speedup graphs for each test case
    for test_case in test_cases:
        # Remove duplicates from the worker list while preserving order
        unique_workers = []
        for w in worker_counts[test_case]:
            if w not in unique_workers:
                unique_workers.append(w)
        
        plot_speedup(results, test_case, unique_workers)
    
    # Display detailed results in terminal
    for test_case in test_cases:
        if results[test_case]["sequential"] is None:
            continue
            
        print(f"\nDetailed Results for {test_case}:")
        print(f"{'Workers':<10} {'Time (s)':<15} {'Speedup':<10}")
        print("-" * 35)
        print(f"{'Sequential':<10} {results[test_case]['sequential']:<15.6f} {1.0:<10.2f}")
        
        # Display results in the original order specified
        for w in worker_counts[test_case]:
            # Skip duplicates
            if w in results[test_case]["parallel"] and results[test_case]["parallel"][w] is not None:
                speedup = results[test_case]["sequential"] / results[test_case]["parallel"][w]
                print(f"{w:<10} {results[test_case]['parallel'][w]:<15.6f} {speedup:<10.2f}")
                # Remove this entry to avoid duplicates
                results[test_case]["parallel"].pop(w, None)
    
    # Save the results to a file
    with open("performance_results.txt", "w") as f:
        for test_case in test_cases:
            f.write(f"Test Case: {test_case}\n")
            if results[test_case]["sequential"] is not None:
                f.write(f"Sequential Time: {results[test_case]['sequential']:.6f} seconds\n")
                f.write("Parallel Results:\n")
                f.write(f"{'Workers':<10} {'Time (s)':<15} {'Speedup':<10}\n")
                f.write("-" * 35 + "\n")
                
                # Remove duplicates while preserving order
                tested_workers = []
                for w in worker_counts[test_case]:
                    if w not in tested_workers and w in results[test_case]["parallel"] and results[test_case]["parallel"][w] is not None:
                        tested_workers.append(w)
                        speedup = results[test_case]["sequential"] / results[test_case]["parallel"][w]
                        f.write(f"{w:<10} {results[test_case]['parallel'][w]:<15.6f} {speedup:<10.2f}\n")
            else:
                f.write("Sequential implementation failed, no results available.\n")
            
            f.write("\n")
    
    print("\nAll results saved to performance_results.txt")

if __name__ == "__main__":
    main()
