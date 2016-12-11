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

#define MAX_NOT_IMPROVING_ITERATIONS 10
#define MAX_ITERATIONS_TO_RESTORE 3

using namespace std;

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
    S.usersCell = new int**[nCells];
    B.usersCell = new int**[nCells];
    R.usersCell = new int**[nCells];
    problem.usersCell = new int**[nCells];
    for (int i = 0; i < this->nCells; i++) {
        S.usersCell[i] = new int*[nCells];
        B.usersCell[i] = new int*[nCells];
        R.usersCell[i] = new int*[nCells];
        problem.usersCell[i] = new int*[nCells];
        for (int m = 0; m < this->nCustomerTypes; m++) {
            S.usersCell[i][m] = new int[nCells];
            B.usersCell[i][m] = new int[nCells];
            R.usersCell[i][m] = new int[nCells];
            problem.usersCell[i][m] = new int[nCells];
        }
    }
    this->bestSolution = new double[4];

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
                S.usersCell[i][m][t] = atoi(word.c_str());
                R.usersCell[i][m][t] = atoi(word.c_str());
                B.usersCell[i][m][t] = atoi(word.c_str());
                problem.usersCell[i][m][t] = atoi(word.c_str());
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
    int iterationsWithoutImprovements = 0;
	int iterationsToRestore = 0;
    int iterations = 0;

    clock_t tStart = clock();

    for (int i = 0; i < nCells; i++) {
        if (problem.activities[i] > 0) {
            Cell tmpCell = Cell(i, problem.activities[i], 0);
            S.cells.push_back(tmpCell);
            R.cells.push_back(tmpCell);
            B.cells.push_back(tmpCell);
        }
    }

    S.cellsWithActivitiesLeft = S.cells.size();
    S.objFunc = 0;
    B.cellsWithActivitiesLeft = B.cells.size();
    B.objFunc = 0;
    R.cellsWithActivitiesLeft = R.cells.size();
    R.objFunc = 0;

    solveGreedy(); //some initial solution

    cout << "INITIAL SOLUTION: " << R.objFunc << "\n";

    copyDataStructure(&B, &R);

    copyDataStructure(&S, &R);

    int k = 0;
    double T = 500;
    double T0 = T;
    double alpha = 0.99985;
    int suboptimalMovements = 0;
    vector<double> Ttrace;
    Ttrace.push_back(T);

    double time = ((clock() - tStart) / (double) CLOCKS_PER_SEC );

    while (time < 4.92) {
        
		copyDataStructure(&R, &S);

        gentlemanAgreement();

        T *= alpha;
        Ttrace.push_back(T);

        if (R.objFunc < S.objFunc) {
            copyDataStructure(&S, &R);
		}

		if (iterationsWithoutImprovements == MAX_NOT_IMPROVING_ITERATIONS) {
            iterationsToRestore++;
            if (iterationsToRestore == MAX_ITERATIONS_TO_RESTORE) {
                copyDataStructure(&S, &B);
                iterationsToRestore = 0;
                iterationsWithoutImprovements = 0;
            } else {
                gentlemenClub(floor(R.cells.size() * 0.4));
                iterationsWithoutImprovements = 0;
                if (R.objFunc < S.objFunc) {
                    copyDataStructure(&S, &R);
                } else if (R.objFunc > S.objFunc) {
                    double p = exp(-((R.objFunc - S.objFunc)) / T);
                    initRandSeed;
                    double r = (((double) (rand() % 9)) / 10.0) + 0.1;
                    if (r < p) {
                        copyDataStructure(&S, &R);
                    }
                }
            }
        }

        if (S.objFunc < B.objFunc) {
            copyDataStructure(&B, &S);
            iterationsWithoutImprovements = 0;
			iterationsToRestore = 0;
        }  else iterationsWithoutImprovements++;

        time = ((clock() - tStart) / (double) CLOCKS_PER_SEC );

        iterations++;
    }

    cout << "ITERATIONS: " << iterations << "\n";


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

void Heuristic::solveGreedy() {

    for (int x = 0; x < R.cells.size(); x++) {
		int demand = R.cells[x].activities;
        bool tooPrecise = true;
		bool continueSearching = true;
        while (demand > 0 && continueSearching) {
            bool agentFound = false;
            double bestRatio = INT32_MAX;
            int best_J, best_M, best_T;
            for (int j = 0; j < nCells; j++) {
                for (int m = 0; m < nCustomerTypes; m++) {
                    for (int t = 0; t < nTimeSteps; t++) {
						if ((R.cells[x].i != j) && (R.usersCell[j][m][t] != 0)) {
							if (tooPrecise) {
								if (R.cells[x].activities >= problem.n[m]) {
									double ratio = problem.costs[j][R.cells[x].i][m][t] / problem.n[m];
									if (ratio < bestRatio) {
										agentFound = true;
										bestRatio = ratio;
										best_J = j;
										best_M = m;
										best_T = t;
									}
								}
							} else {
								double ratio = problem.costs[j][R.cells[x].i][m][t] / demand;
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
                int i = R.cells[x].i;
                int j = best_J;
                int m = best_M;
                int t = best_T;

                if (demand > (problem.n[m] *  R.usersCell[j][m][t])) {
                    R.cells[x].partialObjFunc += R.usersCell[j][m][t] * problem.costs[j][i][m][t];
                    R.objFunc += R.usersCell[j][m][t] * problem.costs[j][i][m][t];
                    demand -= R.usersCell[j][m][t] * problem.n[m];
                    Agent tmpAgent = Agent(j, m, t, R.usersCell[j][m][t]);
                    R.cells[x].usedAgents.push_back(tmpAgent);
                    R.usersCell[j][m][t] = 0;
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
                    R.cells[x].partialObjFunc += howmany * problem.costs[j][i][m][t];
                    R.objFunc += howmany * problem.costs[j][i][m][t];
                    demand -= howmany * problem.n[m];
                    R.cells[x].usedAgents.push_back(tmpAgent);
                    R.usersCell[j][m][t] -= howmany;
                }
            } else if (!agentFound && tooPrecise) {
				tooPrecise = false;
			} else if (!agentFound && !tooPrecise){
				continueSearching = false;
			}
            R.cells[x].activities = demand;
        }
    }

}

void Heuristic::gentlemanGreedy() {

    for (int x = 0; x < R.cells.size(); x++) {
        int demand = R.cells[x].activities;
        bool continueSearching = true;
        while (demand > 0 && continueSearching) {
            bool agentFound = false;
            double bestRatio = INT32_MAX;
            int best_J, best_M, best_T;
            for (int j = 0; j < nCells && (demand > 0); j++) {
                for (int m = 0; m < nCustomerTypes && (demand > 0); m++) {
                    for (int t = 0; t < nTimeSteps && (demand > 0); t++) {
                        if ((R.cells[x].i != j) && (R.usersCell[j][m][t] > 0)) {
                            double ratio;
                            if (demand >= problem.n[m]) {
                                ratio = problem.costs[j][R.cells[x].i][m][t] / problem.n[m];
                            } else ratio = problem.costs[j][R.cells[x].i][m][t] / demand;
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

				int i = R.cells[x].i;
				int j = best_J;
				int m = best_M;
				int t = best_T;

				if (demand > (problem.n[m] * R.usersCell[j][m][t])) {
					R.cells[x].partialObjFunc += R.usersCell[j][m][t] * problem.costs[j][i][m][t];
					R.objFunc += R.usersCell[j][m][t] * problem.costs[j][i][m][t];
					demand -= R.usersCell[j][m][t] * problem.n[m];
					Agent tmpAgent = Agent(j, m, t, R.usersCell[j][m][t]);
					R.cells[x].usedAgents.push_back(tmpAgent);
					R.usersCell[j][m][t] = 0;
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
					R.cells[x].partialObjFunc += howmany * problem.costs[j][i][m][t];
					R.objFunc += howmany * problem.costs[j][i][m][t];
					demand -= howmany * problem.n[m];
					R.cells[x].usedAgents.push_back(tmpAgent);
					R.usersCell[j][m][t] -= howmany;
				}
			} else {
				continueSearching = false;
            }

            R.cells[x].activities = demand;
        }
    }
}

void Heuristic::gentlemanAgreement() {
	initRandSeed;
	int y;
	int x = rand() % R.cells.size();
	do {
		y = rand() % R.cells.size();
	} while (x == y);

	int beneficiary, donor;

	if (R.cells[x].partialObjFunc > R.cells[y].partialObjFunc) {
		donor = y;
		beneficiary = x;
	}
	else {
		donor = x;
		beneficiary = y;
	}

	int i_d = R.cells[donor].i;
	int j_d = R.cells[donor].usedAgents.front().j;
	int m_d = R.cells[donor].usedAgents.front().m;
	int t_d = R.cells[donor].usedAgents.front().t;
	int n_d = R.cells[donor].usedAgents.front().n;

	int i_b = R.cells[beneficiary].i;
	int j_b = R.cells[beneficiary].usedAgents.front().j;
	int m_b = R.cells[beneficiary].usedAgents.front().m;
	int t_b = R.cells[beneficiary].usedAgents.front().t;
	int n_b = R.cells[beneficiary].usedAgents.front().n;

	R.objFunc -= problem.costs[j_d][i_d][m_d][t_d] * n_d;
	R.cells[donor].partialObjFunc -= problem.costs[j_d][i_d][m_d][t_d] * n_d;
	R.usersCell[j_d][m_d][t_d] += n_d;
	R.cells[beneficiary].usedAgents.front().n = 0;

	R.objFunc -= problem.costs[j_b][i_b][m_b][t_b] * n_b;
	R.cells[beneficiary].partialObjFunc -= problem.costs[j_b][i_b][m_b][t_b] * n_b;
	R.usersCell[j_b][m_b][t_b] += n_b;

	R.cells[beneficiary].activities += problem.n[m_b] * n_b;
	R.cells[beneficiary].usedAgents.pop_front();

	gentlemanGreedy();

    R.cells[donor].activities += problem.n[m_d] * n_d;
    R.cells[donor].usedAgents.pop_front();

    gentlemanGreedy();

    if (R.cells[beneficiary].activities > 0 || R.cells[donor].activities > 0) {
        copyDataStructure(&R, &S);
    }

}

void Heuristic::gentlemenClub(int n) {
    vector<int> chosenOnes, activitiesToAdd;

    initRandSeed;

    for (int i = 0; i < n; i++) {
        int x;
        do {
            x = rand() % R.cells.size();
        } while (find(chosenOnes.begin(), chosenOnes.end(), x) != chosenOnes.end());
        activitiesToAdd.push_back(0);
        chosenOnes.push_back(x);
    }


    for (int l = 0; l < n; l++) {
        int iterations = R.cells[chosenOnes[l]].usedAgents.size();

        for (int k = 0; k < iterations; k++) {
            int i = R.cells[chosenOnes[l]].i;
            int j = R.cells[chosenOnes[l]].usedAgents.front().j;
            int m = R.cells[chosenOnes[l]].usedAgents.front().m;
            int t = R.cells[chosenOnes[l]].usedAgents.front().t;
            int howmany = R.cells[chosenOnes[l]].usedAgents.front().n;

            R.objFunc -= problem.costs[j][i][m][t] * howmany;
            R.cells[chosenOnes[l]].partialObjFunc -= problem.costs[j][i][m][t] * howmany;
            R.usersCell[j][m][t] += howmany;
            activitiesToAdd[l] += problem.n[m] * howmany;
            R.cells[chosenOnes[l]].usedAgents.pop_front();
        }
    }

    for (int i = 0; i < n; i++) {
        R.cells[chosenOnes[i]].activities += activitiesToAdd[i];
        gentlemanGreedy();
		if (R.cells[chosenOnes[i]].activities > 0) {
			copyDataStructure(&R, &S);
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
