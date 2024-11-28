#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#define MAX_VILLES 100
#define MAX_ELEPHANTS 100

typedef struct {
    double x, y;
} Point;

double distance(Point a, Point b) {
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

// Function for Clan Updating Operator(elephant position)
void clanUpdatingOperator(double* newElephantPositions, double* oldElephantPositions, double matriarchPosition, double a, int numElephants) {
    for (int j = 0; j < numElephants; ++j) {
        double r = ((double)rand() / RAND_MAX);
        double newPosition = oldElephantPositions[j] + a * (matriarchPosition - oldElephantPositions[j]) * r;
        newElephantPositions[j] = newPosition;
    }
}

// Function to update the matriarch position using Clan Updating Operator
double updateMatriarchPosition(double* elephantPositions, int numElephants, double b) {
    double sum = 0.0;
    for (int j = 0; j < numElephants; ++j) {
        sum += elephantPositions[j];
    }
    return b * (sum / numElephants);
}

// Function to calculate the fitness of an elephant (example implementation)
double calculateFitness(double position) {
    return position * position; //  fitness is the square of the position
}

// Implement the Separating Operator
void separatingOperator(double* elephantPositions, int numElephants, double xmin, double xmax) {
    for (int i = 0; i < numElephants; ++i) {
        // Identify the elephant with the worst fitness in each clan
        int worstIndex = i; // Assume the current elephant is the worst for now
        double worstFitness = calculateFitness(elephantPositions[worstIndex]);
        for (int j = i + 1; j < numElephants; ++j) {
            double currentFitness = calculateFitness(elephantPositions[j]);
            if (currentFitness < worstFitness) {
                worstFitness = currentFitness;
                worstIndex = j;
            }
        }

        // Replace the elephant with the worst fitness with a new randomly generated position
        double r = ((double)rand() / RAND_MAX);
        elephantPositions[worstIndex] = xmin + (xmax - xmin) * r;
    }
}

// Function to calculate the total distance for a given order of cities
double calculateTotalDistance(int n, Point villes[MAX_VILLES], int* order) {
    double totalDistance = 0.0;

    for (int i = 0; i < n - 1; ++i) {
        totalDistance += distance(villes[order[i]], villes[order[i + 1]]);
    }

    // Add distance from the last city back to the starting city
    totalDistance += distance(villes[order[n - 1]], villes[order[0]]);

    return totalDistance;
}

// Function to display the order of cities
void displayCityOrder(int n, int* order, double totalDistance) {
    printf("Chemin optimal : Ville %d", order[0] + 1);
    for (int i = 1; i < n; ++i) {
        printf(" -> Ville %d", order[i] + 1);
    }
    printf(" -> Ville %d\n", order[0] + 1);
    printf("Distance totale : %.2f\n", totalDistance);
}

// Function to perform 2-opt local search
void twoOpt(int n, Point villes[MAX_VILLES], int* order) {
    double bestDistance = calculateTotalDistance(n, villes, order);
    int improvement = 1;

    while (improvement) {
        improvement = 0;
        for (int i = 1; i < n - 1; ++i) {
            for (int j = i + 1; j < n; ++j) {
                int* newOrder = (int*)malloc(n * sizeof(int));
                for (int k = 0; k < i; ++k) {
                    newOrder[k] = order[k];
                }
                for (int k = i, l = j; k < j + 1; ++k, --l) {
                    newOrder[k] = order[l];
                }
                for (int k = j + 1; k < n; ++k) {
                    newOrder[k] = order[k];
                }

                double newDistance = calculateTotalDistance(n, villes, newOrder);
                if (newDistance < bestDistance) {
                    bestDistance = newDistance;
                    for (int k = 0; k < n; ++k) {
                        order[k] = newOrder[k];
                    }
                    improvement = 1;
                }
                free(newOrder);
            }
        }
    }
}

void elephantHerdingOptimization(int n, Point villes[MAX_VILLES], int* order) {
    const double a = 0.1; // Facteur de mise à jour du clan
    const double b = 0.1; // Facteur de mise à jour de la matriarche
    const int numElephants = MAX_ELEPHANTS;
    const int MaxGen = 100; // Nombre maximum de générations
    const int nClan = 5; // Nombre de clans (vous pouvez ajuster cela en fonction de votre problème)

    double elephantPositions[MAX_ELEPHANTS];
    double matriarchPosition = ((double)rand() / RAND_MAX);

    for (int i = 0; i < numElephants; ++i) {
        elephantPositions[i] = ((double)rand() / RAND_MAX);
    }

    double bestDistance = INFINITY;
    int bestOrder[MAX_VILLES];

    int t = 0; // Compteur de génération
    while (t < MaxGen) {
        matriarchPosition = updateMatriarchPosition(elephantPositions, numElephants, b);

        // Mise à jour des positions des éléphants par opérateur de mise à jour du clan
        for (int ci = 0; ci < nClan; ++ci) {
            clanUpdatingOperator(elephantPositions, elephantPositions, matriarchPosition, a, numElephants);
        }

        int currentOrder[MAX_VILLES];
        for (int i = 0; i < n; ++i) {
            currentOrder[i] = i;
        }

        twoOpt(n, villes, currentOrder);

        double currentDistance = calculateTotalDistance(n, villes, currentOrder);



        if (currentDistance < bestDistance) {
            bestDistance = currentDistance;
            for (int i = 0; i < n; ++i) {
                bestOrder[i] = currentOrder[i];
            }
        }

        // Opérateur de séparation
        for (int ci = 0; ci < nClan; ++ci) {
            separatingOperator(elephantPositions, numElephants, 0, 1);
        }

        t++;
    }

    for (int i = 0; i < n; ++i) {
        order[i] = bestOrder[i];
    }
}





int main() {
    srand((unsigned int)time(NULL));

    int n;
    Point villes[MAX_VILLES];

    // Read Berlin52 instance from file
    FILE* file = fopen("berlin52.tsp", "r");
    if (file == NULL) {
        printf("Error opening file.\n");
        return 1;
    }

    char line[100];
    while (fgets(line, sizeof(line), file) != NULL) {
        if (strstr(line, "NODE_COORD_SECTION") != NULL) {
            break;
        }
    }

    int id;
    double x, y;
    int i = 0;
    while (fscanf(file, "%d %lf %lf", &id, &x, &y) == 3) {
        villes[i].x = x;
        villes[i].y = y;
        i++;
    }
    n = i;

    fclose(file);

    int order[MAX_VILLES];
    elephantHerdingOptimization(n, villes, order);
    double totalDistance = calculateTotalDistance(n, villes, order);

    printf("Chemin optimal trouve par l'algorithme de l'optimisation de l'herd d'éléphants avec 2-opt :\n");
    displayCityOrder(n, order, totalDistance);

    return 0;
}
