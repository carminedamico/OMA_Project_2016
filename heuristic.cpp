#include <iostream>
#include <random>
#include "heuristic.h"
#include <time.h>
#include <algorithm>

#ifdef _WIN32
#include <Windows.h>
	#include <algorithm>
	#define initRandSeed SYSTEMTIME st; GetSystemTime(&st); srand(st.wMilliseconds);
#else
#include <ctime>
#include <unistd.h>
#define initRandSeed struct timespec spec; clock_gettime(CLOCK_REALTIME,&spec); srand(spec.tv_nsec);
#endif

using namespace std;

int MAX_NOT_IMPROVING_ITERATIONS = 0;
int MAX_ITERATIONS_TO_RESTORE = 0;

int HIVE_SIZE;
int EMPLOYED_BEES;
int ONLOOKER_BEES_PER_EMPLOYED;
int FREE_BEES;

struct less_than_key {
	inline bool operator() (const Bee& struct1, const Bee& struct2) {
		return (struct1.solution.objFunc < struct2.solution.objFunc);
	}
};

Heuristic::Heuristic(string path){
    this->hasSolution = false;
    string line;
    string word;

    ifstream iffN(path.c_str());

    if (!iffN.is_open()) {
        cout << "Impossible to open" << path << endl;
        cin.get();
        exit(1);
    }

    getline(iffN, line);
    std::replace(line.begin(), line.end(), ';', ' ');
    istringstream iss(line);
    iss >> word;
    this->nCells = atoi(word.c_str());
    iss >> word;
    this->nTimeSteps = atoi(word.c_str());
    iss >> word;
    this->nCustomerTypes = atoi(word.c_str());

	if (nCells > 200) {
		MAX_NOT_IMPROVING_ITERATIONS = 100;
		HIVE_SIZE = 13;
		EMPLOYED_BEES = 3;
		ONLOOKER_BEES_PER_EMPLOYED = 3;
		MAX_ITERATIONS_TO_RESTORE = 3;
	}
	else if (nCells > 50) {
		MAX_NOT_IMPROVING_ITERATIONS = 20;
		HIVE_SIZE = 100;
		EMPLOYED_BEES = 4;
		ONLOOKER_BEES_PER_EMPLOYED = 10;
		MAX_ITERATIONS_TO_RESTORE = 3;
	}
	else {
		MAX_NOT_IMPROVING_ITERATIONS = 3;
		HIVE_SIZE = 100;
		EMPLOYED_BEES = 5;
		ONLOOKER_BEES_PER_EMPLOYED = 10;
		MAX_ITERATIONS_TO_RESTORE = 5;
	}

    // Memory allocation
    solution = new int***[nCells];
    problem.costs = new double***[nCells];
    for (int i = 0; i < this->nCells; i++) {
        problem.costs[i] = new double**[nCells];
        solution[i] = new int**[nCells];
        for (int j = 0; j < this->nCells; j++) {
            problem.costs[i][j] = new double*[nCustomerTypes];
            solution[i][j] = new int*[nCustomerTypes];
            for (int m = 0; m < this->nCustomerTypes; m++) {
                problem.costs[i][j][m] = new double[nTimeSteps];
                solution[i][j][m] = new int[nTimeSteps];
            }
        }
    }
    problem.n = new int[nCustomerTypes];
    problem.activities = new int[nCells];
    B.usersCell = new int**[nCells];
    problem.usersCell = new int**[nCells];
    for (int i = 0; i < this->nCells; i++) {
        B.usersCell[i] = new int*[nCustomerTypes];
        problem.usersCell[i] = new int*[nCustomerTypes];
        for (int m = 0; m < this->nCustomerTypes; m++) {
            B.usersCell[i][m] = new int[nTimeSteps];
            problem.usersCell[i][m] = new int[nTimeSteps];
        }
    }
    this->bestSolution = new double[4];

	for (int i = 0; i < HIVE_SIZE; i++) {
		Bee employed;
		/*for (int j = 0; j < ONLOOKER_BEES_PER_EMPLOYED; j++) {
			Bee onlooker;
			employed.onlookers.push_back(onlooker);
		}*/
		bees.push_back(employed);
	}

	for (int b = 0; b < bees.size(); b++) {
		bees[b].solution.usersCell = new int**[nCells];
		for (int i = 0; i < this->nCells; i++) {
			bees[b].solution.usersCell[i] = new int*[nCells];
			for (int m = 0; m < this->nCustomerTypes; m++) {
				bees[b].solution.usersCell[i][m] = new int[nCells];
			}
		}
	}

    getline(iffN, line);
    getline(iffN, line);
    std::replace(line.begin(), line.end(), ';', ' ');
    istringstream issN(line);
    for (int m = 0; m < nCustomerTypes; m++) {
        issN >> word;
        problem.n[m] = atoi(word.c_str());
    }

    getline(iffN, line);
    for (int m = 0; m < nCustomerTypes; m++) {
        for (int t = 0; t < nTimeSteps; t++) {
            getline(iffN, line);// linea con m e t
            for (int i = 0; i < nCells; i++) {
                getline(iffN, line);// linea della matrice c_{ij} per t ed m fissati
                istringstream issC(line);
                for (int j = 0; j < nCells; j++) {
                    issC >> word;
                    problem.costs[i][j][m][t] = atoi(word.c_str());
                }
            }
        }
    }

    getline(iffN, line);
    getline(iffN, line);
    std::replace(line.begin(), line.end(), ';', ' ');
    istringstream issA(line);
    for (int i = 0; i < nCells; i++) {
        issA >> word;
        problem.activities[i] = atoi(word.c_str());
    }

    getline(iffN, line);
    for (int m = 0; m < nCustomerTypes; m++) {
        for (int t = 0; t < nTimeSteps; t++) {
            getline(iffN, line);
            getline(iffN, line);
            std::replace(line.begin(), line.end(), ';', ' ');
            istringstream issU(line);
            for (int i = 0; i < nCells; i++) {
                issU >> word;
                B.usersCell[i][m][t] = atoi(word.c_str());
                problem.usersCell[i][m][t] = atoi(word.c_str());
            }
        }
    }

	for (int m = 0; m < nCustomerTypes; m++) {
		for (int t = 0; t < nTimeSteps; t++) {
			for (int i = 0; i < nCells; i++) {
				for (int b = 0; b < bees.size(); b++) {
					bees[b].solution.usersCell[i][m][t] = B.usersCell[i][m][t];
				}
			}
		}
	}

    string best_path = "./Optimal_Solutions.csv";

    ifstream iffB(best_path.c_str());

    if (!iffB.is_open()) {
        cout << "Impossible to open" << "./Optimal_Solutions.csv" << endl;
        cin.get();
        bestSolutionKnown = false;
    } else bestSolutionKnown = true;

    if (bestSolutionKnown) {
        do {
            getline(iffB, line);
            std::replace(line.begin(), line.end(), ';', ' ');
            istringstream isb(line);
            isb >> word;
            word += '.';
        } while(path.find(word) == string::npos);

        istringstream isb(line);
        isb >> word;
        isb >> word;
        for (int k = 0; k < 4; k++) {
            isb >> word;
            this -> bestSolution[k] = atof(word.c_str());
        }
    }

}

void Heuristic::copyDataStructure(Solution* Y, Solution* X) {
    (*Y).objFunc = (*X).objFunc;
    (*Y).cellsWithActivitiesLeft = (*X).cellsWithActivitiesLeft;

    (*Y).cells.clear();
    vector<Cell>().swap(((*Y).cells));

    for (int i = 0; i < (*X).cells.size(); i++) {
        Cell tmpCell = Cell((*X).cells[i].i, (*X).cells[i].activities,
							(*X).cells[i].partialObjFunc);
        (*Y).cells.push_back(tmpCell);

        for (int j = 0; j < (*X).cells[i].usedAgents.size(); j++) {
            Agent tmpAgent = Agent((*X).cells[i].usedAgents[j].j,
								(*X).cells[i].usedAgents[j].m,
                                (*X).cells[i].usedAgents[j].t,
								(*X).cells[i].usedAgents[j].n);
            (*Y).cells[i].usedAgents.push_back(tmpAgent);
        }

    }

    for (int j = 0; j < nCells; j++)
        for (int m = 0; m < nCustomerTypes; m++)
            for (int t = 0; t < nTimeSteps; t++)
                (*Y).usersCell[j][m][t] = (*X).usersCell[j][m][t];

}

void Heuristic::Metaheuristic(vector<double>& stat) {

    clock_t tStart = clock();	// Time start

	/* Data structures initialization */
	for (int i = 0; i < nCells; i++) {
		if (problem.activities[i] > 0) {
			Cell tmpCell = Cell(i, problem.activities[i], 0);
			B.cells.push_back(tmpCell);
			for (int j = 0; j < HIVE_SIZE; ++j) {
				bees[j].solution.cells.push_back(tmpCell);
				bees[j].solution.objFunc = 0;
			}
		}
	}

	/* Find some initial solution for each bee in the hive */
	for (int i = 0; i < HIVE_SIZE; ++i) {
		solveGreedy(&bees[i].solution);
	}

	/* Order the hive from best bee to worst bee. */
	sort(bees.begin(), bees.end(), less_than_key());

    /*  After the ordering, the first bee in the hive is the best solution actually found */
	copyDataStructure(&B, &bees[0].solution);

    cout << "\tINITIAL SOLUTION: " << B.objFunc << "\n";

	int iterations = 0;	// Initialize iterations counter

    double time = ((clock() - tStart) / (double) CLOCKS_PER_SEC );

    double averageIterationTime = 0;

	while (time < (4.9 - averageIterationTime)) {

		// Ordering of the employed bees (EMPLOYED BEES ONLY!)
		sort(bees.begin(), bees.begin()+EMPLOYED_BEES, less_than_key());

		/* Reclutation of the onlooker bees. */
		int reclutationIndex = HIVE_SIZE - 1;   // Starting from the last bee in the hive (the worst)

		for (int i = 0; i < EMPLOYED_BEES; ++i) {	// For each employed bee (for each of the first EMPLOYED_BEES bees in the hive)
			bees[i].onlookers.clear();				// Remove the previous onlooker bees (to avoid cuncurrency in the reclutation of the onlookers)
			for (int j = reclutationIndex; j > reclutationIndex - ONLOOKER_BEES_PER_EMPLOYED + i; j--) {
				bees[i].onlookers.push_back(&bees[j]);	// Assign an onlooker to the emplouyed bee
			}
			reclutationIndex -= (ONLOOKER_BEES_PER_EMPLOYED + i);
		}

		/* For each employed bee */
		for (int b = 0; b < EMPLOYED_BEES; b++) {

			/* Each onlooker bee perform an exploration in the neighborhood of the employed bee she is assigned. The index of the best onlooker
				bee is saved and used to verify the quality of the exlporation.*/
			int bestOnlookerIndex = 0;	// The index of the best onlooker bee after the exploration

			for (int o = 0; o < bees[b].onlookers.size(); o++) {						// For each onlooker
				copyDataStructure(&bees[b].onlookers[o]->solution, &bees[b].solution);	// The onlooker bee starts from the position of the employed bee
				gentlemenClub(&bees[b].onlookers[o]->solution, &bees[b].solution, 2);	// The onlooker performs a gentlemanClub
				if (bees[b].onlookers[o]->solution.objFunc < bees[b].onlookers[bestOnlookerIndex]->solution.objFunc) {
					bestOnlookerIndex = o;		// Eventually save the index of the onlooker if she has the best obj function
				}
			}

            for (int o = 0; o < bees[b].onlookers.size(); o++) {						// For each onlooker
                gentlemanAgreement(&bees[b].onlookers[o]->solution, &bees[b].solution);	// The onlooker performs a gentlemanAgreement
                if (bees[b].onlookers[o]->solution.objFunc < bees[b].onlookers[bestOnlookerIndex]->solution.objFunc) {
                    bestOnlookerIndex = o;		// Eventually save the index of the onlooker if she has the best obj function
                }
            }

			/* Is the objective function of the best onlooker better than the objective function of the employed? */
			if (bees[b].onlookers[bestOnlookerIndex]->solution.objFunc < bees[b].solution.objFunc) {
				copyDataStructure(&bees[b].solution, &bees[b].onlookers[bestOnlookerIndex]->solution);	// If yes, move the employed in the new solution

			} else {
				bees[b].iterationsWithoutImprovements++;	// otherwise increment the iterationWithoutImprovements counter
			}

			/* If in the last MAX_NOT_IMPROVING_ITERATIONS iterations there was no improvement in the employed bee objective function,
				maybe the employed bee is in a local optimum. To be sure to move in a good point of the space (where she is sure there will be
				enought food), she reaches the most promising zone: the zone where there is the best solution.
				This movement is done only if the employed bee actually is not in the best solution zone. */
			if (bees[b].iterationsWithoutImprovements == MAX_NOT_IMPROVING_ITERATIONS) {
				if (bees[b].solution.objFunc > B.objFunc) {		// If the employed bee is not in the best solution zone
					copyDataStructure(&bees[b].solution, &B);	// she moves to the best solution
				}
				else {											// Otherwise
					bees[b].iterationsToRestore++;				// she increment the iterationToRestoreCounter
					/* If iteration to restore counter is equal to MAX_ITERATION_TO_RESTORE, it maybe the employed bee explored enought
						the zone; to avoid the possibility remain in a local optimum, she move away to search another solution.
						If the employed bee is making the wrong choice, she will return to the best solution in the future.
						Furthermore other employed bees will remain in the beest solution zone to explore it, so the choice to move
						away is not critical. */
					if (bees[b].iterationsToRestore == MAX_ITERATIONS_TO_RESTORE) {
						gentlemenClub(&bees[b].solution, &B, floor(bees[b].solution.cells.size()*0.9));	// The movement is performed with a 60% gentlemanClub
						bees[b].iterationsToRestore = 0;											    // After the movement, reset iterationToRestore
					}
				}
				bees[b].iterationsWithoutImprovements = 0;	// Reset iterationWithoutImprovements
			}

			/* Verify if the employed bee found a new best solution. */
			if (bees[b].solution.objFunc < B.objFunc) {
				copyDataStructure(&B, &bees[b].solution);
				/* If the bee has found a new best solution, maybe she is on the right way. Reset the counters to avoid an undesired movement in the
					upcoming future */
				bees[b].iterationsWithoutImprovements = 0;
				bees[b].iterationsToRestore = 0;
			}

			time = ((clock() - tStart) / (double)CLOCKS_PER_SEC);

			iterations++;

            averageIterationTime = time / iterations;

		}
    }

    cout << "\tITERATIONS: " << iterations << "\n";


    stat.push_back(B.objFunc);
    stat.push_back(((clock() - tStart) / (double) CLOCKS_PER_SEC ));

    hasSolution=true;

    for (int i = 0; i < nCells; i++)
        for (int j = 0; j < nCells; j++)
            for (int m = 0; m < nCustomerTypes; m++)
                for (int t = 0; t < nTimeSteps; t++)
                    solution[i][j][m][t] = 0;

    for (int n = 0; n < B.cells.size(); n++) {
        for (int l = 0; l < B.cells[n].usedAgents.size(); l++) {
            int i = B.cells[n].i;
            int j = B.cells[n].usedAgents[l].j;
            int m = B.cells[n].usedAgents[l].m;
            int t = B.cells[n].usedAgents[l].t;
            int used = B.cells[n].usedAgents[l].n;
            solution[j][i][m][t] += used;
        }
    }

}

void Heuristic::solveGreedy(Solution *R) {

	random_shuffle(R->cells.begin(), R->cells.end());

    for (int x = 0; x < R->cells.size(); x++) {
		int demand = R->cells[x].activities;
        bool tooPrecise = true;
		bool continueSearching = true;
        while (demand > 0 && continueSearching) {
            bool agentFound = false;
            double bestRatio = INT32_MAX;
            int best_J, best_M, best_T;
            for (int j = 0; j < nCells; j++) {
                for (int m = 0; m < nCustomerTypes; m++) {
                    for (int t = 0; t < nTimeSteps; t++) {
						if ((R->cells[x].i != j) && (R->usersCell[j][m][t] != 0)) {
							if (tooPrecise) {
								if (R->cells[x].activities >= problem.n[m]) {
									double ratio = problem.costs[j][R->cells[x].i][m][t] / problem.n[m];
									if (ratio < bestRatio) {
										agentFound = true;
										bestRatio = ratio;
										best_J = j;
										best_M = m;
										best_T = t;
									}
								}
							} else {
								double ratio = problem.costs[j][R->cells[x].i][m][t] / demand;
								if (ratio < bestRatio) {
									agentFound = true;
									bestRatio = ratio;
									best_J = j;
									best_M = m;
									best_T = t;
								}
							}
                        }
                    }
                }
            }

            if (agentFound) {
                int i = R->cells[x].i;
                int j = best_J;
                int m = best_M;
                int t = best_T;

                if (demand > (problem.n[m] *  R->usersCell[j][m][t])) {
                    R->cells[x].partialObjFunc += R->usersCell[j][m][t] * problem.costs[j][i][m][t];
                    R->objFunc += R->usersCell[j][m][t] * problem.costs[j][i][m][t];
                    demand -= R->usersCell[j][m][t] * problem.n[m];
                    Agent tmpAgent = Agent(j, m, t, R->usersCell[j][m][t]);
                    R->cells[x].usedAgents.push_back(tmpAgent);
                    R->usersCell[j][m][t] = 0;
                } else {
                    int howmany;
                    if (tooPrecise) {
                        howmany = floor(demand / problem.n[m]);
                    } else {
                        if ( demand % problem.n[m] == 0) {
                            howmany = floor(demand / problem.n[m]);
                        } else {
                            howmany = floor(demand / problem.n[m]) + 1;
                        }
                    }
                    Agent tmpAgent = Agent(j, m, t, howmany);
                    R->cells[x].partialObjFunc += howmany * problem.costs[j][i][m][t];
                    R->objFunc += howmany * problem.costs[j][i][m][t];
                    demand -= howmany * problem.n[m];
                    R->cells[x].usedAgents.push_back(tmpAgent);
                    R->usersCell[j][m][t] -= howmany;
                }
            } else if (!agentFound && tooPrecise) {
				tooPrecise = false;
			} else if (!agentFound && !tooPrecise){
				continueSearching = false;
			}
            R->cells[x].activities = demand;
        }
    }

}

void Heuristic::gentlemanGreedy(Solution *R) {

    for (int x = 0; x < R->cells.size(); x++) {
        int demand = R->cells[x].activities;
        bool continueSearching = true;
        while (demand > 0 && continueSearching) {
            bool agentFound = false;
            double bestRatio = INT32_MAX;
            int best_J, best_M, best_T;
            for (int j = 0; j < nCells && (demand > 0); j++) {
                for (int m = 0; m < nCustomerTypes && (demand > 0); m++) {
                    for (int t = 0; t < nTimeSteps && (demand > 0); t++) {
                        if ((R->cells[x].i != j) && (R->usersCell[j][m][t] > 0)) {
                            double ratio;
                            if (demand >= problem.n[m]) {
                                ratio = problem.costs[j][R->cells[x].i][m][t] / problem.n[m];
                            } else ratio = problem.costs[j][R->cells[x].i][m][t] / demand;
                            if (ratio < bestRatio) {
                                agentFound = true;
                                bestRatio = ratio;
                                best_J = j;
                                best_M = m;
                                best_T = t;
                            }
                        }
                    }
                }
            }

			if (agentFound) {

				int i = R->cells[x].i;
				int j = best_J;
				int m = best_M;
				int t = best_T;

				if (demand > (problem.n[m] * R->usersCell[j][m][t])) {
					R->cells[x].partialObjFunc += R->usersCell[j][m][t] * problem.costs[j][i][m][t];
					R->objFunc += R->usersCell[j][m][t] * problem.costs[j][i][m][t];
					demand -= R->usersCell[j][m][t] * problem.n[m];
					Agent tmpAgent = Agent(j, m, t, R->usersCell[j][m][t]);
					R->cells[x].usedAgents.push_back(tmpAgent);
					R->usersCell[j][m][t] = 0;
				}
				else {
					int howmany;
					if (demand % problem.n[m] == 0) {
						howmany = floor(demand / problem.n[m]);
					}
					else {
						howmany = floor(demand / problem.n[m]) + 1;
					}
					Agent tmpAgent = Agent(j, m, t, howmany);
					R->cells[x].partialObjFunc += howmany * problem.costs[j][i][m][t];
					R->objFunc += howmany * problem.costs[j][i][m][t];
					demand -= howmany * problem.n[m];
					R->cells[x].usedAgents.push_back(tmpAgent);
					R->usersCell[j][m][t] -= howmany;
				}
			} else {
				continueSearching = false;
            }

            R->cells[x].activities = demand;
        }
    }
}

void Heuristic::gentlemanAgreement(Solution *R, Solution *S) {
	initRandSeed;
	int y;
	int x = rand() % R->cells.size();
	do {
		y = rand() % R->cells.size();
	} while (x == y);

	int beneficiary, donor;

	if (R->cells[x].partialObjFunc > R->cells[y].partialObjFunc) {
		donor = y;
		beneficiary = x;
	}
	else {
		donor = x;
		beneficiary = y;
	}

	int i_d = R->cells[donor].i;
	int j_d = R->cells[donor].usedAgents.front().j;
	int m_d = R->cells[donor].usedAgents.front().m;
	int t_d = R->cells[donor].usedAgents.front().t;
	int n_d = R->cells[donor].usedAgents.front().n;

	int i_b = R->cells[beneficiary].i;
	int j_b = R->cells[beneficiary].usedAgents.front().j;
	int m_b = R->cells[beneficiary].usedAgents.front().m;
	int t_b = R->cells[beneficiary].usedAgents.front().t;
	int n_b = R->cells[beneficiary].usedAgents.front().n;

	R->objFunc -= problem.costs[j_d][i_d][m_d][t_d] * n_d;
	R->cells[donor].partialObjFunc -= problem.costs[j_d][i_d][m_d][t_d] * n_d;
	R->usersCell[j_d][m_d][t_d] += n_d;
	R->cells[beneficiary].usedAgents.front().n = 0;

	R->objFunc -= problem.costs[j_b][i_b][m_b][t_b] * n_b;
	R->cells[beneficiary].partialObjFunc -= problem.costs[j_b][i_b][m_b][t_b] * n_b;
	R->usersCell[j_b][m_b][t_b] += n_b;

	R->cells[beneficiary].activities += problem.n[m_b] * n_b;
	R->cells[beneficiary].usedAgents.pop_front();

	gentlemanGreedy(R);

    R->cells[donor].activities += problem.n[m_d] * n_d;
    R->cells[donor].usedAgents.pop_front();

    gentlemanGreedy(R);

    if (R->cells[beneficiary].activities > 0 || R->cells[donor].activities > 0) {
        copyDataStructure(R, S);
    }

}

void Heuristic::gentlemenClub(Solution *R, Solution *S, int n) {
    vector<int> chosenOnes, activitiesToAdd, randomNumbers;

    for (int i = 0; i < R->cells.size(); ++i) {
        randomNumbers.push_back(i);
    }

    random_shuffle(randomNumbers.begin(), randomNumbers.end());

    for (int i = 0; i < n; i++) {
        activitiesToAdd.push_back(0);
        chosenOnes.push_back(randomNumbers[i]);
    }


    for (int l = 0; l < n; l++) {
        int iterations = R->cells[chosenOnes[l]].usedAgents.size();

        for (int k = 0; k < iterations; k++) {
            int i = R->cells[chosenOnes[l]].i;
            int j = R->cells[chosenOnes[l]].usedAgents.front().j;
            int m = R->cells[chosenOnes[l]].usedAgents.front().m;
            int t = R->cells[chosenOnes[l]].usedAgents.front().t;
            int howmany = R->cells[chosenOnes[l]].usedAgents.front().n;

            R->objFunc -= problem.costs[j][i][m][t] * howmany;
            R->cells[chosenOnes[l]].partialObjFunc -= problem.costs[j][i][m][t] * howmany;
            R->usersCell[j][m][t] += howmany;
            activitiesToAdd[l] += problem.n[m] * howmany;
            R->cells[chosenOnes[l]].usedAgents.pop_front();
        }
    }

    for (int i = 0; i < n; i++) {
        R->cells[chosenOnes[i]].activities += activitiesToAdd[i];
        gentlemanGreedy(R);
		if (R->cells[chosenOnes[i]].activities > 0) {
			copyDataStructure(R, S);
			break;
		}
    }

}

void Heuristic::writeKPI(string path, string instanceName, vector<double> stat){
    if (!hasSolution)
        return;

    ofstream fileO(path, ios::app);
    if(!fileO.is_open())
        return;

    fileO << instanceName << ";" << stat[0] << ";" << stat[1];
    for(int i=2; i<stat.size(); i++)
        fileO <<  ";" << stat[i];

	fileO << ";" << (this->bestSolution[0]);
    fileO << endl;

    fileO.close();

}

void Heuristic::writeSolution(string path) {
    if (!hasSolution)
        return;

    ofstream fileO(path);
    if(!fileO.is_open())
        return;

    fileO << this->nCells << "; " << this->nTimeSteps << "; " << this->nCustomerTypes << endl;
    for (int m = 0; m < this->nCustomerTypes; m++)
        for (int t = 0; t < this->nTimeSteps; t++)
            for (int i = 0; i < this->nCells; i++)
                for (int j = 0; j < this->nCells; j++)
                    if (solution[i][j][m][t] > 0)
                        fileO << i << ";" << j << ";" << m << ";" << t << ";"
							<< solution[i][j][m][t] << endl;

    fileO.close();
}

eFeasibleState Heuristic::isFeasible(string path) {

    string line;
    string word;
    int nCellsN;
    int nTimeStepsN;
    int nCustomerTypesN;
    int i, j, m, t;


    ifstream iffN(path.c_str());

    if (!iffN.is_open()) {
        cout << "Impossible to open" << path << endl;
        exit(1);
    }

    getline(iffN, line);
    std::replace(line.begin(), line.end(), ';', ' ');
    istringstream iss(line);
    iss >> word; // nCells
    nCellsN = atoi(word.c_str());
    iss >> word; // nTimeSteps
    nTimeStepsN = atoi(word.c_str());
    iss >> word; // nCustomerTypes
    nCustomerTypesN = atoi(word.c_str());

    int**** solutionN = new int***[nCells];
    for (i = 0; i < nCellsN; i++) {
        solutionN[i] = new int**[nCells];
        for (j = 0; j < nCellsN; j++) {
            solutionN[i][j] = new int*[nCustomerTypes];
            for (m = 0; m < nCustomerTypesN; m++) {
                solutionN[i][j][m] = new int[nTimeSteps];
                for ( t = 0; t < nTimeStepsN; t++) {
                    solutionN[i][j][m][t] = 0;
                }
            }
        }
    }

    while (getline(iffN, line)) {
        std::replace(line.begin(), line.end(), ';', ' ');
        istringstream iss(line);
        iss >> word; // i
        i = atoi(word.c_str());
        iss >> word; // j
        j = atoi(word.c_str());
        iss >> word; // m
        m = atoi(word.c_str());
        iss >> word; // t
        t = atoi(word.c_str());
        iss >> word; // value
        solutionN[i][j][m][t] = atoi(word.c_str());
    }

    // Demand
    bool feasible = true;
    int expr;
    for (int i = 0; i < nCells; i++) {
        expr = 0;
        for (int j = 0; j < nCells; j++)
            for (int m = 0; m < nCustomerTypes; m++)
                for (int t = 0; t < nTimeSteps; t++)
                    expr += problem.n[m] * solutionN[j][i][m][t];
        if (expr < problem.activities[i])

            feasible = false;
    }

    if (!feasible)
        return NOT_FEASIBLE_DEMAND;

    // Max Number of users
    for (int i = 0; i < nCells; i++) {
        for (int m = 0; m < nCustomerTypes; m++) {
            for (int t = 0; t < nTimeSteps; t++) {
                expr = 0;
                for (int j = 0; j < nCells; j++)
                    expr += solutionN[i][j][m][t];
                if (expr > problem.usersCell[i][m][t]) {
                    feasible = false;
                }
            }
        }
    }

    if(!feasible)
        return NOT_FEASIBLE_USERS;

    return FEASIBLE;
}

void Heuristic::getStatSolution(vector<double>& stat) {
    if (!hasSolution)
        return;

    int* tipi = new int[nCustomerTypes];
    for (int m = 0; m < nCustomerTypes; m++)
        tipi[m] = 0;

    for (int i = 0; i < nCells; i++)
        for (int j = 0; j < nCells; j++)
            for (int t = 0; t < nTimeSteps; t++)
                for (int m = 0; m < nCustomerTypes; m++)
                    if (solution[i][j][m][t] > 0)
                        tipi[m] += solution[i][j][m][t];
    for (int m = 0; m < nCustomerTypes; m++)
        stat.push_back(tipi[m]);

}
