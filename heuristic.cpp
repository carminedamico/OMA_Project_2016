//
// Created by Luca Gobbato on 04/10/16.
//

#include <iostream>
#include <random>
#include "heuristic.h"
#include <time.h>
#include <ctime>
#include <unistd.h>

#ifdef _WIN32
#include <Windows.h>
	#include <algorithm>
	#define initRandSeed SYSTEMTIME st; GetSystemTime(&st); srand(st.wMilliseconds);
#else
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
    problem.costs = new double***[nCells];
    for (int i = 0; i < this->nCells; i++) {
        problem.costs[i] = new double**[nCells];
        for (int j = 0; j < this->nCells; j++) {
            problem.costs[i][j] = new double*[nCustomerTypes];
            for (int m = 0; m < this->nCustomerTypes; m++) {
                problem.costs[i][j][m] = new double[nTimeSteps];
            }
        }
    }
    problem.n = new int[nCustomerTypes];
    problem.activities = new int[nCells];
    problem.usersCell = new int**[nCells];
    for (int i = 0; i < this->nCells; i++) {
        problem.usersCell[i] = new int*[nCustomerTypes];
        for (int m = 0; m < this->nCustomerTypes; m++) {
            problem.usersCell[i][m] = new int[nTimeSteps];
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
            }
        }
    }

    for (int i = 0; i < nCells; i++) {
       if (problem.activities[i] > 0) {
           Cell tmpCell = Cell(i, problem.activities[i]);
           this -> cells.push_back(tmpCell);
       }
   }

    string best_path = "./Optimal_Solutions.csv";

    ifstream iffB(best_path.c_str());

    if (!iffB.is_open()) {
        cout << "Impossible to open" << "./Optimal_Solutions.csv" << endl;
        cin.get();
        exit(1);
    }

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

void Heuristic::Metaheuristic(vector<double>& stat) {
    double sObjFunc = 0, // initial objective function
            bestObjFunc = 0; // best objective function
    clock_t tStart = clock();
    int maxNotImprovingIterations = 30;

    solveGreedy(&sObjFunc, &cells); //some initial solution

    bestObjFunc = sObjFunc;

    vector<Cell> S(cells), R;

    double time = ((clock() - tStart) / (double) CLOCKS_PER_SEC );

    int k = 0;
    double T = 500;
    double T0 = T;
    double alpha = 0.9995;
    int suboptimalMovements = 0;
    vector<double> Ttrace;

    while (time < 4.9) {
        if (++k % 1000 == 0) {
            Ttrace.push_back(T);
            //T = T0*pow(alpha,k);
            T -= T*alpha;
        }

        R.clear();
        R = S;
        double rObjFunc = sObjFunc;
        gentlemanAgreement(&rObjFunc, &R);

        if (rObjFunc < sObjFunc) {
            S.clear();
            S = R;
            sObjFunc = rObjFunc;
        } else if (rObjFunc > sObjFunc) {
            double p = exp(-(rObjFunc-sObjFunc)/T);
            initRandSeed;
            double r = ((double)(rand() % 10)) / 10.0;

            if (r < p) {
                S = R;
                sObjFunc = rObjFunc;
                suboptimalMovements++;
            }

        }

        if (sObjFunc < bestObjFunc) {
            cells.clear();
            cells = S;
            bestObjFunc = sObjFunc;
            maxNotImprovingIterations = 30;
        } else {
            //maxNotImprovingIterations --;
            //if (maxNotImprovingIterations == 0) {
            //    gentlemenClub(&sObjFunc, &S, (int) cells.size() / 2);
            //    maxNotImprovingIterations = 30;
            //}
        }

        time = ((clock() - tStart) / (double) CLOCKS_PER_SEC );

    }

    stat.push_back(bestObjFunc);
    stat.push_back(((clock() - tStart) / (double) CLOCKS_PER_SEC ));

    hasSolution=true;

}

void Heuristic::solveGreedy(double* ObjFunc, vector<Cell>* cells) {

    for (int x = 0; x < (*cells).size(); x++) {
        bool AgentFound = true;
        bool tooPrecise = true;
        while ((*cells)[x].activities > 0) {
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

                if ((*cells)[x].activities > (problem.n[m] *  problem.usersCell[j][m][t])) {
                    (*cells)[x].partialObjFunc += problem.usersCell[j][m][t] * problem.costs[j][i][m][t];
                    *ObjFunc += problem.usersCell[j][m][t] * problem.costs[j][i][m][t];
                    (*cells)[x].activities -= problem.usersCell[j][m][t] * problem.n[m];
                    Agent tmpAgent = Agent(j, m, t, problem.usersCell[j][m][t]);
                    (*cells)[x].usedAgents.push_back(tmpAgent);
                    problem.usersCell[j][m][t] = 0;
                } else {
                    int howmany = floor((*cells)[x].activities / problem.n[m]);
                    if (tooPrecise) {
                        howmany = floor((*cells)[x].activities / problem.n[m]);
                    } else {
                        if ( demand % problem.n[m] == 0) {
                            howmany = floor((*cells)[x].activities / problem.n[m]);
                        } else {
                            howmany = floor((*cells)[x].activities / problem.n[m]) + 1;
                        }
                    }
                    Agent tmpAgent = Agent(j, m, t, howmany);
                    (*cells)[x].partialObjFunc += howmany * problem.costs[j][i][m][t];
                    *ObjFunc += howmany * problem.costs[j][i][m][t];
                    (*cells)[x].activities -= howmany * problem.n[m];
                    (*cells)[x].usedAgents.push_back(tmpAgent);
                    problem.usersCell[j][m][t] -= howmany;
                }
            }

            (*cells)[x].activities = demand;
        }
    }
}

void Heuristic::gentlemanGreedy(double *ObjFunc, vector<Cell> *cells) {

    for (int x = 0; x < (*cells).size(); x++) {
        while ((*cells)[x].activities > 0) {
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

            if ((*cells)[x].activities > (problem.n[m] *  problem.usersCell[j][m][t])) {

                (*cells)[x].partialObjFunc += problem.usersCell[j][m][t] * problem.costs[j][i][m][t];
                *ObjFunc += problem.usersCell[j][m][t] * problem.costs[j][i][m][t];
                (*cells)[x].activities -= problem.usersCell[j][m][t] * problem.n[m];
                Agent tmpAgent = Agent(j, m, t, problem.usersCell[j][m][t]);
                (*cells)[x].usedAgents.push_back(tmpAgent);
                problem.usersCell[j][m][t] = 0;

            } else {

                int howmany;
                if ( demand % problem.n[m] == 0) {
                    howmany = floor((*cells)[x].activities / problem.n[m]);
                } else {
                    howmany = floor((*cells)[x].activities / problem.n[m]) + 1;
                }
                Agent tmpAgent = Agent(j, m, t, howmany);
                (*cells)[x].partialObjFunc += howmany * problem.costs[j][i][m][t];
                *ObjFunc += howmany * problem.costs[j][i][m][t];
                (*cells)[x].activities -= howmany * problem.n[m];
                (*cells)[x].usedAgents.push_back(tmpAgent);
                problem.usersCell[j][m][t] -= howmany;

            }
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

    if ((*cells)[x].partialObjFunc > (*cells)[x].partialObjFunc) {
        donor = y;
        beneficiary = x;
    } else {
        donor = x;
        beneficiary = y;
    }

    int i_d = (*cells)[donor].i;
    int j_d = (*cells)[donor].usedAgents[0].j;
    int m_d = (*cells)[donor].usedAgents[0].m;
    int t_d = (*cells)[donor].usedAgents[0].t;
    int n_d = (*cells)[donor].usedAgents[0].n;

    int i_b = (*cells)[beneficiary].i;
    int j_b = (*cells)[beneficiary].usedAgents[0].j;
    int m_b = (*cells)[beneficiary].usedAgents[0].m;
    int t_b = (*cells)[beneficiary].usedAgents[0].t;
    int n_b = (*cells)[beneficiary].usedAgents[0].n;

    *ObjFunc -= problem.costs[j_d][i_d][m_d][t_d] * n_d;
    (*cells)[donor].partialObjFunc -= problem.costs[j_d][i_d][m_d][t_d] * n_d;
    problem.usersCell[j_d][m_d][t_d] += n_d;
    (*cells)[donor].usedAgents.pop_front();

    *ObjFunc -= problem.costs[j_b][i_b][m_b][t_b] * n_b;
    (*cells)[beneficiary].partialObjFunc -= problem.costs[j_b][i_b][m_b][t_b] * n_b;
    problem.usersCell[j_b][m_b][t_b] += n_d;
    (*cells)[beneficiary].activities += problem.n[m_b] * n_b;
    (*cells)[beneficiary].usedAgents.pop_front();

    gentlemanGreedy(ObjFunc, cells);

    usleep(500);

    (*cells)[donor].activities += problem.n[(*cells)[donor].usedAgents[0].m];

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
    for (int i = 0; i < this->cells.size(); i++) {
        for (int j = 0; j < cells[i].usedAgents.size(); j++) {
            fileO << cells[i].usedAgents[j].j << ";" << cells[i].i << ";" << cells[i].usedAgents[j].m << ";" << cells[i].usedAgents[j].t endl;
        }
    }

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
                if (expr > problem.usersCell[i][m][t])
                    feasible = false;
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

    for (int i = 0; i < cells.size(); i++) {
        for (int j = 0; j < cells[i].usedAgents.size(); j++) {
            tipi[cells[i].usedAgents[j].m] += cells[i].usedAgents[j].n;
        }
    }
    for (int m = 0; m < nCustomerTypes; m++)
        stat.push_back(tipi[m]);

}
