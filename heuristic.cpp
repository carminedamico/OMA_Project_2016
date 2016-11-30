//
// Created by Luca Gobbato on 04/10/16.
//

#include <iostream>
#include <random>
#include "heuristic.h"
#include <unistd.h>

#ifdef _WIN32
#include <Windows.h>
	#include <algorithm>
	#define initRandSeed SYSTEMTIME st; GetSystemTime(&st); srand(st.wMilliseconds);
#else
#include <time.h>
#include <ctime>
#define initRandSeed struct timespec spec; clock_gettime(CLOCK_REALTIME,&spec); srand(spec.tv_nsec);
#endif

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
    problem.usersCell = new int**[nCells];
    problem.ORIGINALusersCell = new int**[nCells];
    for (int i = 0; i < this->nCells; i++) {
        problem.usersCell[i] = new int*[nCustomerTypes];
        problem.ORIGINALusersCell[i] = new int*[nCustomerTypes];
        for (int m = 0; m < this->nCustomerTypes; m++) {
            problem.usersCell[i][m] = new int[nTimeSteps];
            problem.ORIGINALusersCell[i][m] = new int[nTimeSteps];
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
                problem.usersCell[i][m][t] = atoi(word.c_str());
                problem.ORIGINALusersCell[i][m][t] = problem.usersCell[i][m][t];
            }
        }
    }

    for (int i = 0; i < nCells; i++) {
       if (problem.activities[i] > 0) {
           Cell tmpCell = Cell(i, problem.activities[i], 0);
           this -> cells.push_back(tmpCell);
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

void Heuristic::copyDataStructure(vector<Cell>* Y, vector<Cell>* X) {
    (*Y).clear();
    vector<Cell>().swap((*Y));

    for (int i = 0; i < (*X).size(); i++) {
        Cell tmpCell = Cell((*X)[i].i, (*X)[i].activities, (*X)[i].partialObjFunc);
        (*Y).push_back(tmpCell);

        int iterations = (*X)[i].usedAgents.size();
        for (int j = 0; j < iterations; j++) {
            Agent tmpAgent = Agent((*X)[i].usedAgents.front().j, (*X)[i].usedAgents.front().m,
                                   (*X)[i].usedAgents.front().t, (*X)[i].usedAgents.front().n);
            (*X)[i].usedAgents.pop();
            (*X)[i].usedAgents.push(tmpAgent);
            (*Y)[i].usedAgents.push(tmpAgent);
        }
    }

}

void Heuristic::Metaheuristic(vector<double>& stat) {
    double sObjFunc = 0, // initial objective function
         bestObjFunc = 0, // best objective function
            rObjFunc = 0; // tmp
    clock_t tStart = clock();
    int maxNotImprovingIterations = 30;

    solveGreedy(&bestObjFunc, &cells); //some initial solution

    sObjFunc = bestObjFunc;

    rObjFunc = bestObjFunc;

    cout << bestObjFunc << "\n";

    vector<Cell> S;

    copyDataStructure(&S, &cells);


    for (int i1 = 0; i1 < 5000; i1++) {
        vector<Cell> R;
        copyDataStructure(&R, &S);
        gentlemanAgreement(&rObjFunc, &R);
    }

    double time = ((clock() - tStart) / (double) CLOCKS_PER_SEC );

    int k = 0;
    double T = 500;
    double T0 = T;
    double alpha = 0.9995;
    int suboptimalMovements = 0;
    vector<double> Ttrace;

    /*while (time < 4.9) {
        if (++k % 1000 == 0) {
            Ttrace.push_back(T);
            T -= T*alpha;
        }


        copyDataStructure(&S, &R);

        gentlemanAgreement(&sObjFunc, &S);

        if (sObjFunc < rObjFunc) {
            copyDataStructure(&R, &S);
            rObjFunc = sObjFunc;
        } else if (sObjFunc > rObjFunc) {
            double p = exp(-(rObjFunc-sObjFunc)/T);
            initRandSeed;
            double r = ((double)(rand() % 10)) / 10.0;

            if (r < p) {
                copyDataStructure(&R, &S);
                rObjFunc = sObjFunc;
            }

        }

        if (rObjFunc < bestObjFunc) {
            copyDataStructure(&cells, &R);
            bestObjFunc = rObjFunc;
            //maxNotImprovingIterations = 30;
        } else {
            //maxNotImprovingIterations --;
            //if (maxNotImprovingIterations == 0) {
            //    gentlemenClub(&sObjFunc, &S, (int) cells.size() / 2);
            //    maxNotImprovingIterations = 30;
            //}
        }

        time = ((clock() - tStart) / (double) CLOCKS_PER_SEC );

    }*/

    stat.push_back(bestObjFunc);
    stat.push_back(((clock() - tStart) / (double) CLOCKS_PER_SEC ));

    hasSolution=true;

    for (int i = 0; i < nCells; i++)
        for (int j = 0; j < nCells; j++)
            for (int m = 0; m < nCustomerTypes; m++)
                for (int t = 0; t < nTimeSteps; t++)
                    solution[i][j][m][t] = 0;

    for (int n = 0; n < cells.size(); n++) {
        int iterations = cells[n].usedAgents.size();
        for (int l = 0; l < iterations; l++) {
            int i = cells[n].i;
            int j = cells[n].usedAgents.front().j;
            int m = cells[n].usedAgents.front().m;
            int t = cells[n].usedAgents.front().t;
            int used = cells[n].usedAgents.front().n;
            solution[j][i][m][t] += used;
            cells[n].usedAgents.pop();
        }

    }

}

void Heuristic::solveGreedy(double* ObjFunc, vector<Cell>* cells) {

    for (int x = 0; x < (*cells).size(); x++) {
        int demand = (*cells)[x].activities;
        bool AgentFound = true;
        bool tooPrecise = true;
        while (demand > 0) {
            if (!AgentFound) tooPrecise = false;
            AgentFound = false;
            double bestRatio = 9999999;
            int bestj, bestm, bestt;
            for (int j = 0; j < nCells && (demand > 0); j++) {
                for (int m = 0; m < nCustomerTypes && (demand > 0); m++) {
                    if (tooPrecise) {
                        if (demand < problem.n[m]) {
                            break;
                        }
                    }
                    for (int t = 0; t < nTimeSteps && (demand > 0); t++) {
                        if (((*cells)[x].i != j) && (problem.usersCell[j][m][t] != 0)) {
                            double ratio = problem.costs[j][(*cells)[x].i][m][t] / problem.n[m];
                            if (ratio < bestRatio) {
                                AgentFound = true;
                                bestRatio = ratio;
                                bestj = j;
                                bestm = m;
                                bestt = t;
                            }
                        }
                    }
                }
            }

            if (AgentFound) {
                int i = (*cells)[x].i;
                int j = bestj;
                int m = bestm;
                int t = bestt;

                if (demand > (problem.n[m] *  problem.usersCell[j][m][t])) {
                    (*cells)[x].partialObjFunc += problem.usersCell[j][m][t] * problem.costs[j][i][m][t];
                    *ObjFunc += problem.usersCell[j][m][t] * problem.costs[j][i][m][t];
                    demand -= problem.usersCell[j][m][t] * problem.n[m];
                    Agent tmpAgent = Agent(j, m, t, problem.usersCell[j][m][t]);
                    (*cells)[x].usedAgents.push(tmpAgent);
                    problem.usersCell[j][m][t] = 0;
                } else {
                    int howmany = floor(demand / problem.n[m]);
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
                    (*cells)[x].partialObjFunc += howmany * problem.costs[j][i][m][t];
                    *ObjFunc += howmany * problem.costs[j][i][m][t];
                    demand -= howmany * problem.n[m];
                    (*cells)[x].usedAgents.push(tmpAgent);
                    problem.usersCell[j][m][t] -= howmany;
                }
            }

            (*cells)[x].activities = demand;
        }
    }
}

void Heuristic::gentlemanGreedy(double *ObjFunc, vector<Cell> *cells) {

    for (int x = 0; x < (*cells).size(); x++) {
        int demand = (*cells)[x].activities;
        while (demand > 0) {
            double bestRatio = 9999999;
            int bestj, bestm, bestt;
            for (int j = 0; j < nCells && (demand > 0); j++) {
                for (int m = 0; m < nCustomerTypes && (demand > 0); m++) {
                    for (int t = 0; t < nTimeSteps && (demand > 0); t++) {
                        if (((*cells)[x].i != j) && (problem.usersCell[j][m][t] != 0)) {
                            double ratio = problem.costs[j][(*cells)[x].i][m][t] / problem.n[m];
                            if (ratio < bestRatio) {
                                bestRatio = ratio;
                                bestj = j;
                                bestm = m;
                                bestt = t;
                            }
                        }
                    }
                }
            }

            int i = (*cells)[x].i;
            int j = bestj;
            int m = bestm;
            int t = bestt;

            if (demand > (problem.n[m] *  problem.usersCell[j][m][t])) {
                (*cells)[x].partialObjFunc += problem.usersCell[j][m][t] * problem.costs[j][i][m][t];
                *ObjFunc += problem.usersCell[j][m][t] * problem.costs[j][i][m][t];
                demand -= problem.usersCell[j][m][t] * problem.n[m];
                Agent tmpAgent = Agent(j, m, t, problem.usersCell[j][m][t]);
                (*cells)[x].usedAgents.push(tmpAgent);
                problem.usersCell[j][m][t] = 0;
            } else {
                int howmany;
                if ( demand % problem.n[m] == 0 ) {
                    howmany = floor(demand / problem.n[m]);
                } else {
                    howmany = floor(demand / problem.n[m]) + 1;
                }
                Agent tmpAgent = Agent(j, m, t, howmany);
                (*cells)[x].partialObjFunc += howmany * problem.costs[j][i][m][t];
                *ObjFunc += howmany * problem.costs[j][i][m][t];
                demand -= howmany * problem.n[m];
                (*cells)[x].usedAgents.push(tmpAgent);
                problem.usersCell[j][m][t] -= howmany;
            }
            (*cells)[x].activities = demand;
        }
    }
}

void Heuristic::gentlemanAgreement(double* ObjFunc, vector<Cell>* cells) {
    initRandSeed;

    int x = rand() % (*cells).size();
    int y;
    do {
        y = rand() % (*cells).size();
    } while ( x == y);

    int beneficiary, donor;

    if ((*cells)[x].partialObjFunc > (*cells)[y].partialObjFunc) {
        donor = y;
        beneficiary = x;
    } else {
        donor = x;
        beneficiary = y;
    }

    int i_d = (*cells)[donor].i;
    int j_d = (*cells)[donor].usedAgents.front().j;
    int m_d = (*cells)[donor].usedAgents.front().m;
    int t_d = (*cells)[donor].usedAgents.front().t;
    int n_d = (*cells)[donor].usedAgents.front().n;

    int i_b = (*cells)[beneficiary].i;
    int j_b = (*cells)[beneficiary].usedAgents.front().j;
    int m_b = (*cells)[beneficiary].usedAgents.front().m;
    int t_b = (*cells)[beneficiary].usedAgents.front().t;
    int n_b = (*cells)[beneficiary].usedAgents.front().n;

    *ObjFunc -= problem.costs[j_d][i_d][m_d][t_d] * n_d;
    (*cells)[donor].partialObjFunc -= problem.costs[j_d][i_d][m_d][t_d] * n_d;
    problem.usersCell[j_d][m_d][t_d] += n_d;
    if (problem.usersCell[j_d][m_d][t_d] > problem.ORIGINALusersCell[j_d][m_d][t_d]) {
        cout << "1 " << j_d << " " << problem.usersCell[j_d][m_d][t_d] << " ORIGINALE " << problem.ORIGINALusersCell[j_d][m_d][t_d] << "\n";
    }
    (*cells)[beneficiary].usedAgents.front().n = 0;

    *ObjFunc -= problem.costs[j_b][i_b][m_b][t_b] * n_b;
    (*cells)[beneficiary].partialObjFunc -= problem.costs[j_b][i_b][m_b][t_b] * n_b;
    problem.usersCell[j_b][m_b][t_b] += n_b;
    if (problem.usersCell[j_b][m_b][t_b] > problem.ORIGINALusersCell[j_b][m_b][t_b]) {
        cout << "2 " << j_b << " " << problem.usersCell[j_b][m_b][t_b] << " ORIGINALE " << problem.ORIGINALusersCell[j_b][m_b][t_b] << "\n";
    }
    (*cells)[beneficiary].activities += problem.n[m_b] * n_b;
    (*cells)[beneficiary].usedAgents.pop();

    gentlemanGreedy(ObjFunc, cells);

    usleep(500);

    (*cells)[donor].activities += problem.n[m_d] * n_d;
    (*cells)[donor].usedAgents.pop();

    gentlemanGreedy(ObjFunc, cells);

}

/*void Heuristic::gentlemenClub(double* ObjFunc, vector<Cell>* cells, int n) {
    vector<int> chosenOnes, activitiesToAdd;

    initRandSeed;

    for (int i = 0; i < n; i++) {
        int x;
        do {
            x = rand() % (*cells).size();
        } while (find(chosenOnes.begin(), chosenOnes.end(), x) != chosenOnes.end());
        activitiesToAdd.push_back(0);
        chosenOnes.push_back(x);
    }

    //order by crescent objFun?

    for (int i = 0; i < n; i++) {
        while ((*cells)[chosenOnes[i]].usedAgents.size() > 0) {
            *ObjFunc -= problem.costs[(*cells)[chosenOnes[i]].usedAgents[0].j][(*cells)[chosenOnes[i]].i][(*cells)[chosenOnes[i]].usedAgents[0].m][(*cells)[chosenOnes[i]].usedAgents[0].t] * (*cells)[chosenOnes[i]].usedAgents[0].n;
            (*cells)[chosenOnes[i]].partialObjFunc -= problem.costs[(*cells)[chosenOnes[i]].usedAgents[0].j][(*cells)[chosenOnes[i]].i][(*cells)[chosenOnes[i]].usedAgents[0].m][(*cells)[chosenOnes[i]].usedAgents[0].t] * (*cells)[chosenOnes[i]].usedAgents[0].n;
            problem.usersCell[(*cells)[chosenOnes[i]].usedAgents[0].j][(*cells)[chosenOnes[i]].usedAgents[0].m][(*cells)[chosenOnes[i]].usedAgents[0].t] += (*cells)[chosenOnes[i]].usedAgents[0].n;
            activitiesToAdd[i] += problem.n[(*cells)[chosenOnes[i]].usedAgents[0].m] * (*cells)[chosenOnes[i]].usedAgents[0].n;
            (*cells)[chosenOnes[i]].usedAgents[0].n = 0;
            (*cells)[chosenOnes[i]].usedAgents.pop_front();
        }
    }

    for (int i = 0; i < n; i++) {
        (*cells)[chosenOnes[i]].activities += activitiesToAdd[i];
        gentlemanGreedy(ObjFunc, cells);
        usleep(500);
    }

}
*/

void Heuristic::writeKPI(string path, string instanceName, vector<double> stat){
    if (!hasSolution)
        return;

    ofstream fileO(path, ios::app);
    if(!fileO.is_open())
        return;

    fileO << instanceName << ";" << stat[0] << ";" << stat[1];
    for(int i=2; i<stat.size(); i++)
        fileO <<  ";" << stat[i];

    fileO << (((stat[0] - this->bestSolution[0]) /  this->bestSolution[0]) * 100) << "%";
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
                        fileO << i << ";" << j << ";" << m << ";" << t << ";" << solution[i][j][m][t] << endl;

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
    for (int i = 0; i < nCells; i++)
        for (int m = 0; m < nCustomerTypes; m++)
            for (int t = 0; t < nTimeSteps; t++) {
                expr = 0;
                for (int j = 0; j < nCells; j++)
                    expr += solutionN[i][j][m][t];
                if (expr > problem.usersCell[i][m][t]) {
                    cout << i << "\n";
                    feasible = false;
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