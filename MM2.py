import random
import heapq
from collections import deque
import tkinter as tk
from tkinter import ttk  # For themed widgets (optional, looks nicer) pip install -r requirements.txt
from tkinter import messagebox
import math # Needed for isnan check

# --- Core Simulation Logic (modified to be a function) ---

def run_simulation(num_tellers, mean_interarrival, mean_service, sim_duration, seed):
    """Runs the M/M/2 queue simulation and returns results."""

    # --- Parameter Validation (Basic) ---
    if not (num_tellers > 0 and mean_interarrival >= 0 and mean_service >= 0 and sim_duration > 0):
         return {"error": "Invalid parameters. Tellers/Duration must be > 0, Times >= 0."}
    # Allow 0 interarrival time (all arrive at once) but handle division by zero later
    # Allow 0 service time (instant service) but handle division by zero later

    # Use provided seed for reproducibility
    random.seed(seed)

    # --- Simulation State Variables ---
    event_list = []
    sim_time = 0.0
    teller_free_time = [0.0] * num_tellers
    customer_queue = deque()

    # --- Statistics Accumulators ---
    total_wait_time = 0.0
    customers_served_count = 0
    max_queue_length = 0
    teller_total_busy_time = [0.0] * num_tellers
    teller_busy_start_time = {} # {teller_id: start_time}

    # --- Helper Functions (Internal to run_simulation) ---
    def schedule_event(time, event_type, data=None):
        heapq.heappush(event_list, (time, event_type, data))

    def get_random_exponential(mean_time):
        if mean_time <= 0:
            # Return 0 for zero mean, representing instantaneous event
            # Avoid division by zero if mean_time is 0
            return 0.0
        rate = 1.0 / mean_time
        return random.expovariate(rate)

    # --- Event Handler Functions (Internal to run_simulation) ---
    # Using nonlocal to modify variables in the outer function's scope
    def handle_arrival(arrival_time):
        nonlocal max_queue_length, customers_served_count, total_wait_time

        # Schedule next arrival only if mean interarrival > 0
        if mean_interarrival > 0:
             next_arrival_delay = get_random_exponential(mean_interarrival)
             next_arrival_time = arrival_time + next_arrival_delay
             if next_arrival_time < sim_duration:
                 schedule_event(next_arrival_time, "ARRIVAL")
        # If mean_interarrival is 0, all customers arrive effectively at time 0

        idle_teller_id = -1
        for i in range(num_tellers):
            if teller_free_time[i] <= arrival_time:
                idle_teller_id = i
                break

        if idle_teller_id != -1:
            # Serve immediately
            wait_time = 0.0
            total_wait_time += wait_time
            # Increment served count here? No, only when service *completes*
            # The original logic incremented here, let's change to completion for clarity

            service_duration = get_random_exponential(mean_service)
            completion_time = arrival_time + service_duration
            teller_free_time[idle_teller_id] = completion_time
            teller_busy_start_time[idle_teller_id] = arrival_time
            schedule_event(completion_time, "COMPLETION", idle_teller_id)
        else:
            # Queue
            customer_queue.append(arrival_time)
            max_queue_length = max(max_queue_length, len(customer_queue))

    def handle_completion(completion_time, teller_id):
        nonlocal total_wait_time, customers_served_count

        # Finalize busy time calculation for the completed service
        start_time = teller_busy_start_time.pop(teller_id, None) # Remove entry
        if start_time is not None:
             # Ensure start_time makes sense relative to completion_time
             if start_time > completion_time:
                 print(f"Warning: Start time {start_time:.2f} > Completion time {completion_time:.2f} for teller {teller_id}. Clamping busy duration.")
                 busy_duration = 0
             else:
                 busy_duration = completion_time - start_time
             teller_total_busy_time[teller_id] += busy_duration
        else:
            # This might happen if simulation starts with a completion event somehow? Should be rare.
             print(f"Warning: Completion event for Teller {teller_id} at {completion_time:.2f} without a recorded start time.")

        # Increment count *after* service is successfully completed
        customers_served_count += 1

        if customer_queue:
            customer_arrival_time = customer_queue.popleft()
            wait_time = completion_time - customer_arrival_time
            if wait_time < 0:
                # Should not happen, indicates arrival_time > completion_time somehow
                 print(f"Warning: Negative wait time ({wait_time:.2f}) calculated at {completion_time:.2f}. Clamping to 0.")
                 wait_time = 0.0

            total_wait_time += wait_time

            service_duration = get_random_exponential(mean_service)
            next_completion_time = completion_time + service_duration
            teller_free_time[teller_id] = next_completion_time
            teller_busy_start_time[teller_id] = completion_time # Start new busy period now
            schedule_event(next_completion_time, "COMPLETION", teller_id)
        else:
            # Teller becomes idle, free time already set to completion_time
            # No entry in teller_busy_start_time indicates idleness
            pass


    # --- Simulation Setup ---
    # Schedule first arrival if relevant
    if mean_interarrival > 0:
        first_arrival_time = get_random_exponential(mean_interarrival)
        if first_arrival_time < sim_duration:
             schedule_event(first_arrival_time, "ARRIVAL")
    # If mean_interarrival is 0, we might need a different setup logic
    # (e.g., schedule N arrivals at time 0). For simplicity,
    # this MVP assumes mean_interarrival > 0 for standard operation.
    # Handle case of zero service time (completion happens instantly) implicitly.

    # --- Simulation Execution Loop ---
    while event_list:
        try:
            event_time, event_type, data = heapq.heappop(event_list)
        except IndexError:
            break

        if event_time > sim_duration:
            # Put event back? Not needed as we break and use sim_duration as end time.
            break

        sim_time = event_time # Advance clock

        if event_type == "ARRIVAL":
            handle_arrival(event_time)
        elif event_type == "COMPLETION":
            handle_completion(event_time, data)

    # --- Post-Simulation Analysis ---
    final_sim_time = sim_duration

    # Account for tellers still busy at the end
    for teller_id in range(num_tellers):
        if teller_id in teller_busy_start_time:
            start_time = teller_busy_start_time[teller_id]
            # Clamp start_time if it's somehow > final_sim_time
            start_time = min(start_time, final_sim_time)
            final_busy_segment = final_sim_time - start_time
            teller_total_busy_time[teller_id] += final_busy_segment

    # Calculate final statistics
    average_waiting_time = total_wait_time / customers_served_count if customers_served_count > 0 else 0.0
    # Ensure utilization calculation doesn't divide by zero if sim_duration is 0
    teller_utilization = [(busy_time / final_sim_time) * 100 if final_sim_time > 0 else 0 for busy_time in teller_total_busy_time]
    customers_in_queue_at_end = len(customer_queue)

    # Calculate theoretical traffic intensity rho = lambda / (c * mu)
    arrival_rate = 1.0 / mean_interarrival if mean_interarrival > 0 else float('inf') # Infinite arrival rate if time is 0
    service_rate_per_teller = 1.0 / mean_service if mean_service > 0 else float('inf') # Infinite service rate if time is 0
    system_service_rate = num_tellers * service_rate_per_teller

    # Handle edge cases for rho calculation
    if system_service_rate == float('inf'): # Infinite service rate means rho is 0 unless arrival rate is also infinite
        traffic_intensity = 0.0 if arrival_rate != float('inf') else float('nan') # Undefined if both infinite
    elif system_service_rate > 0:
        traffic_intensity = arrival_rate / system_service_rate
    else: # system_service_rate is 0 (because service_time is infinite or tellers=0)
        traffic_intensity = float('inf') if arrival_rate > 0 else 0.0 # Infinite if customers arrive but can't be served

    # Clamp utilization percentages for display
    teller_utilization_clamped = [max(0.0, min(100.0, util)) for util in teller_utilization]


    # Return results as a dictionary
    return {
        "error": None, # Indicate success
        "customers_served": customers_served_count,
        "average_wait": average_waiting_time,
        "max_queue": max_queue_length,
        "final_queue_length": customers_in_queue_at_end,
        "teller_utilization": teller_utilization_clamped,
        "teller_busy_times": teller_total_busy_time, # Raw busy times
        "final_sim_time": final_sim_time,
        "traffic_intensity": traffic_intensity,
    }

# --- GUI Application Class ---

class SimulationApp:
    def __init__(self, master):
        self.master = master
        master.title("M/M/2 Queue Simulator")
        master.geometry("550x550") # Adjusted size

        # Use themed widgets for a slightly better look
        style = ttk.Style()
        style.theme_use('clam') # Or 'alt', 'default', 'classic'

        # Frame for inputs
        input_frame = ttk.LabelFrame(master, text="Simulation Parameters", padding="10")
        input_frame.pack(pady=10, padx=10, fill="x")

        # Input fields using Grid layout within the frame
        ttk.Label(input_frame, text="Number of Tellers:").grid(row=0, column=0, sticky="w", padx=5, pady=2)
        self.tellers_entry = ttk.Entry(input_frame, width=10)
        self.tellers_entry.grid(row=0, column=1, sticky="w", padx=5, pady=2)
        self.tellers_entry.insert(0, "2") # Default value

        ttk.Label(input_frame, text="Mean Interarrival Time:").grid(row=1, column=0, sticky="w", padx=5, pady=2)
        self.arrival_entry = ttk.Entry(input_frame, width=10)
        self.arrival_entry.grid(row=1, column=1, sticky="w", padx=5, pady=2)
        self.arrival_entry.insert(0, "10.0") # Default value

        ttk.Label(input_frame, text="Mean Service Time:").grid(row=2, column=0, sticky="w", padx=5, pady=2)
        self.service_entry = ttk.Entry(input_frame, width=10)
        self.service_entry.grid(row=2, column=1, sticky="w", padx=5, pady=2)
        self.service_entry.insert(0, "15.0") # Default value

        ttk.Label(input_frame, text="Simulation Duration:").grid(row=3, column=0, sticky="w", padx=5, pady=2)
        self.duration_entry = ttk.Entry(input_frame, width=10)
        self.duration_entry.grid(row=3, column=1, sticky="w", padx=5, pady=2)
        self.duration_entry.insert(0, "1000.0") # Default value

        ttk.Label(input_frame, text="Random Seed:").grid(row=4, column=0, sticky="w", padx=5, pady=2)
        self.seed_entry = ttk.Entry(input_frame, width=10)
        self.seed_entry.grid(row=4, column=1, sticky="w", padx=5, pady=2)
        self.seed_entry.insert(0, "42") # Default value

        # Run Button
        self.run_button = ttk.Button(master, text="Run Simulation", command=self.execute_simulation)
        self.run_button.pack(pady=10)

        # Frame for results
        results_frame = ttk.LabelFrame(master, text="Simulation Results", padding="10")
        results_frame.pack(pady=10, padx=10, fill="both", expand=True)

        # Results Text Area
        self.results_text = tk.Text(results_frame, height=15, width=60, wrap=tk.WORD, state=tk.DISABLED) # Start disabled
        # Add Scrollbar
        scrollbar = ttk.Scrollbar(results_frame, orient="vertical", command=self.results_text.yview)
        self.results_text.configure(yscrollcommand=scrollbar.set)

        scrollbar.pack(side=tk.RIGHT, fill="y")
        self.results_text.pack(side=tk.LEFT, fill="both", expand=True)


    def execute_simulation(self):
        """Gets inputs, runs simulation, and displays results."""
        try:
            # Get and validate inputs
            num_tellers = int(self.tellers_entry.get())
            mean_interarrival = float(self.arrival_entry.get())
            mean_service = float(self.service_entry.get())
            sim_duration = float(self.duration_entry.get())
            seed = int(self.seed_entry.get())

            # Basic validation beyond type conversion (handled in run_simulation now)
            if num_tellers <= 0: raise ValueError("Number of tellers must be positive.")
            if mean_interarrival < 0: raise ValueError("Mean interarrival time cannot be negative.")
            if mean_service < 0: raise ValueError("Mean service time cannot be negative.")
            if sim_duration <= 0: raise ValueError("Simulation duration must be positive.")

        except ValueError as e:
            messagebox.showerror("Input Error", f"Invalid input: {e}")
            return

        # Run the simulation
        results = run_simulation(num_tellers, mean_interarrival, mean_service, sim_duration, seed)

        # Check for errors returned from simulation function
        if results.get("error"):
             messagebox.showerror("Simulation Error", results["error"])
             return

        # Format and display results
        self.display_results(num_tellers, mean_interarrival, mean_service, sim_duration, seed, results)


    def display_results(self, num_tellers, mean_interarrival, mean_service, sim_duration, seed, results):
        """Formats and shows results in the text area."""
        # Enable text area for modification
        self.results_text.config(state=tk.NORMAL)
        self.results_text.delete('1.0', tk.END) # Clear previous results

        # Format the output string
        output = f"--- Simulation Parameters Used ---\n"
        output += f"Tellers: {num_tellers}\n"
        output += f"Mean Interarrival Time: {mean_interarrival:.2f}\n"
        output += f"Mean Service Time: {mean_service:.2f}\n"
        output += f"Simulation Duration: {sim_duration:.2f}\n"
        output += f"Random Seed: {seed}\n\n"

        output += f"--- Simulation Results ---\n"
        output += f"Simulation Ended at Time: {results['final_sim_time']:.2f}\n"
        output += f"Total customers served: {results['customers_served']}\n"

        if results['customers_served'] > 0:
            output += f"Average waiting time in queue: {results['average_wait']:.3f}\n"
        else:
            output += "Average waiting time in queue: N/A (no customers served)\n"

        output += f"Maximum queue length observed: {results['max_queue']}\n"
        output += f"Customers in queue at end: {results['final_queue_length']}\n\n"

        output += "Teller Utilization:\n"
        for i, util in enumerate(results['teller_utilization']):
             busy_time = results['teller_busy_times'][i]
             output += f"  Teller {i + 1}: {util:.2f}% (Busy Time: {busy_time:.2f})\n"

        output += f"\nTheoretical Traffic Intensity (rho): "
        rho = results['traffic_intensity']
        if math.isnan(rho):
            output += "NaN (Undefined)\n"
        elif rho == float('inf'):
            output += "Infinity\n"
        else:
            output += f"{rho:.3f}\n"

        if rho >= 1.0 and not math.isinf(rho):
             output += "Warning: Traffic intensity >= 1.0 implies an unstable system (queue likely grows indefinitely).\n"
        elif math.isinf(rho) and rho > 0: # Check specifically for positive infinity
             output += "Warning: Infinite traffic intensity implies an unstable system.\n"


        # Insert the results and disable editing again
        self.results_text.insert(tk.END, output)
        self.results_text.config(state=tk.DISABLED)


# --- Main Execution ---
if __name__ == "__main__":
    root = tk.Tk()
    app = SimulationApp(root)
    root.mainloop()