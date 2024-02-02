import matplotlib.pyplot as plt

# Read data from the text file
file_path = 'pop.txt'  # Replace with the actual path to your file
iterations = []
population_data = []

with open(file_path, 'r') as file:
    for line in file:
        if line.startswith('Total population after iteration'):
            parts = line.split(':')
            iteration_number = int(parts[0].split()[-1])
            population = int(parts[1].strip())
            iterations.append(iteration_number)
            population_data.append(population)

# Plotting
plt.plot(iterations[::100], population_data[::100], marker='o', linestyle='-')

# Annotate the last data point
last_iteration = iterations[-1]
last_population = population_data[-1]
plt.annotate(f'({last_iteration}, {last_population})',
             xy=(last_iteration, last_population),
             xytext=(last_iteration - 1000, last_population + 200),  # Adjust the text position
             arrowprops=dict(facecolor='black', arrowstyle='->'),
             )

# Add a padding to the top
plt.ylim(0, max(population_data) + 500)

plt.title('Population Growth Over Time')
plt.xlabel('Iterations')
plt.ylabel('Population')
plt.grid(True)
plt.show()


# Save the figure as an image file
plt.savefig('population_growth.png')
