import matplotlib.pyplot as plt

# Define the data
total_processes = [1, 2, 4, 8, 16, 32, 64, 128]
max_speedup = [206.133195 / 206.133195, 206.133195 / 103.491485, 206.133195 / 51.930277,
               206.133195 / 26.050067, 206.133195 / 13.020182, 206.133195 / 6.519567,
               206.133195 / 3.341772, 206.133195 / 9.926397]

rounded_speedup = [round(speedup, 2) for speedup in max_speedup]  # Round speedup values

# Plotting the data with rounded speedup values
for i, txt in enumerate(rounded_speedup):
    plt.annotate(txt, (total_processes[i], max_speedup[i]), textcoords="offset points", xytext=(0, 5), ha='center')

plt.plot(total_processes, max_speedup, marker='o', linestyle='-')
plt.xscale('log', base=2)  # Use log scale for better visualization
plt.xticks(total_processes, total_processes)  # Show actual number of processes on x-axis
plt.title('Speedup vs Total Number of Processes/Threads')
plt.xlabel('Total Number of Processes/Threads')
plt.ylabel('Speedup')
# grid only for y-axis
plt.grid(axis='y')
# add small padding to the top of the plot
plt.ylim(0, max(max_speedup) + 5)

# Save the figure as an image file
plt.savefig('speedup_vs_processes.png')

# Show the plot
plt.show()
