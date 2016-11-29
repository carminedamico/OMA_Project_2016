//
// Created by Luca Gobbato on 04/10/16.
//

#ifndef COIOTE_HEURISTIC_HEURISTIC_H
#define COIOTE_HEURISTIC_HEURISTIC_H

#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <list>
#include <map>

using namespace std;

struct Data {
    // the costs ijkt
    double**** costs;
    //number of activity done by one user of type k
    int* n;
    // Activities to be done in node i during the given time horizon
    int* activities;
    // Number of users of type m in i during time period t (ikt)
    int*** usersCell;
};

enum eFeasibleState {
    FEASIBLE,
    NOT_FEASIBLE_DEMAND,
    NOT_FEASIBLE_FLOW,
    NOT_FEASIBLE_USERS
};

class Agent {
public:
    int j;
    int m;
    int t;
    int n;
    int todo;
    Agent(int j, int m, int t, int n) {
        this->j = j;
        this->m = m;
        this->t = t;
        this->n = n;
    }
};

class Cell {
public:
    int i;
    double partialObjFunc;
    int activities;
    deque<Agent> usedAgents;
    Cell(int i, int n) {
        this->i = i;
        this->activities = n;
        this->partialObjFunc = 0;
    }
};



class Heuristic{
private:
    int nTimeSteps;
    int nCustomerTypes;
    int nCells;

    Data problem;

    bool hasSolution;
    double epsilon;
    vector<Cell> cells;

public:
    /**
     * Default constructor
     * @return Heuristic object
     */
    Heuristic(){};
    double* bestSolution;

    /**
     * Constructor from external file
     * @param path path of the external file cotaining the instance of the problem
     * @return
     */
    Heuristic(string path);

    // stat contiene tempo in posizione 0 e objfunction in posizione 1
    void Metaheuristic(vector<double>& stat);

    void solveGreedy(double *ObjFunc, vector<Cell>* cells);

    void gentlemanGreedy(double *ObjFunc, vector<Cell> *cells);

    void gentlemanAgreement(double* ObjFunc, vector<Cell>* cells);

    void gentlemenClub(double* ObjFunc, vector<Cell>* cells, int n);

    void getStatSolution(vector<double>& stat);

    void printProblem();

    void writeKPI(string path, string nameInstance, vector<double> stat);

    void writeSolution(string path);

    eFeasibleState isFeasible(string path);
};

#endif //COIOTE_HEURISTIC_HEURISTIC_H
