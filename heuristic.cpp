//
// Created by Luca Gobbato on 04/10/16.
//

#include <iostream>
#include <random>
#include "heuristic.h"

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
    } while(path.find(word) == string::npos);

    istringstream isb(line);
    isb >> word;
    isb >> word;
    for (int k = 0; k < 4; k++) {
        isb >> word;
        this -> bestSolution[k] = atof(word.c_str());
    }

}

void Heuristic::solveFast(vector<double>& stat, int timeLimit,  bool verbose) {
    double objFun=0;
    clock_t tStart = clock();

    for (int i = 0; i < nCells; i++)
        for (int j = 0; j < nCells; j++)
            for (int m = 0; m < nCustomerTypes; m++)
                for (int t = 0; t < nTimeSteps; t++)
                    solution[i][j][m][t] = 0;


    for (int i = 0; i < nCells; i++) {
        int domanda = problem.activities[i];

        bool nonSoddisfatta = true;
        for (int j = 0; j < nCells && nonSoddisfatta; j++) {
            for (int m = 0; m < nCustomerTypes && nonSoddisfatta; m++) {
                for (int t = 0; t < nTimeSteps && nonSoddisfatta; t++) {
                    if (i != j) {
                        if (domanda > problem.n[m] * problem.usersCell[j][m][t]) {
                            solution[i][j][m][t] = problem.usersCell[j][m][t];
                            problem.usersCell[j][m][t] -= solution[i][j][m][t];
                        }
                        else {
                            if ((domanda%problem.n[m]) == 0) {
                                solution[i][j][m][t] += floor(domanda / problem.n[m]);
                            } else solution[i][j][m][t] += floor(domanda / problem.n[m]) + 1;
                            problem.usersCell[j][m][t] -= solution[i][j][m][t];
                            nonSoddisfatta = false;
                        }
                        if (solution[i][j][m][t] != 0) {
                            objFun += solution[i][j][m][t] * problem.costs[i][j][m][t];
                        }
                        domanda -= problem.n[m]*solution[i][j][m][t];
                    }
                }
            }
        }
    }


    stat.push_back(objFun);
    stat.push_back(((clock() - tStart) / (double) CLOCKS_PER_SEC ));

    hasSolution=true;

}

void Heuristic::solveGreedy(vector<double>& stat, int timeLimit,  bool verbose) {
    double objFun=0;
    clock_t tStart = clock();

    for (int i = 0; i < nCells; i++)
        for (int j = 0; j < nCells; j++)
            for (int m = 0; m < nCustomerTypes; m++)
                for (int t = 0; t < nTimeSteps; t++)
                    solution[i][j][m][t] = 0;

    for (int i = 0; i < nCells; i++) {
        int domanda = problem.activities[i];
        if (domanda != 0) {
            while (domanda > 0) {
                double bestRatio = 9999999;
                int bestj, bestm, bestt;
                for (int j = 0; j < nCells && (domanda > 0); j++) {
                    for (int m = 0; m < nCustomerTypes && (domanda > 0); m++) {
                        for (int t = 0; t < nTimeSteps && (domanda > 0); t++) {
                            if ((i != j) && (problem.usersCell[j][m][t]!=0)) {
                                double ratio = problem.costs[i][j][m][t] / problem.n[m];
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
                if (domanda > problem.n[bestm] * problem.usersCell[bestj][bestm][bestt]) {
                    solution[i][bestj][bestm][bestt] = problem.usersCell[bestj][bestm][bestt];
                    problem.usersCell[bestj][bestm][bestt] -= solution[i][bestj][bestm][bestt];
                }
                else {
                    if ((domanda%problem.n[bestm]) == 0) {
                        solution[i][bestj][bestm][bestt] += floor(domanda / problem.n[bestm]);
                    } else solution[i][bestj][bestm][bestt] += floor(domanda / problem.n[bestm]) + 1;
                    problem.usersCell[bestj][bestm][bestt] -= solution[i][bestj][bestm][bestt];
                }
                if (solution[i][bestj][bestm][bestt] != 0) {
                    objFun += solution[i][bestj][bestm][bestt] * problem.costs[i][bestj][bestm][bestt];
                }
                domanda -= problem.n[bestm]*solution[i][bestj][bestm][bestt];
            }

        }
    }

    stat.push_back(objFun);
    stat.push_back(((clock() - tStart) / (double) CLOCKS_PER_SEC ));

    hasSolution=true;

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

    for (int i = 0; i < nCells; i++)
        for (int j = 0; j < nCells; j++)
            for (int t = 0; t < nTimeSteps; t++)
                for (int m = 0; m < nCustomerTypes; m++)
                    if (solution[i][j][m][t] > 0)
                        tipi[m] += solution[i][j][m][t];
    for (int m = 0; m < nCustomerTypes; m++)
        stat.push_back(tipi[m]);

}