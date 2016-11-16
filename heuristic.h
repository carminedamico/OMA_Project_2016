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


class Heuristic{
private:
    int nTimeSteps;
    int nCustomerTypes;
    int nCells;

    Data problem;

    bool hasSolution;
    double epsilon;
    int**** solution;
    struct move {
        int i;
        int j;
        int m;
        int t;
        int n;
        int* agents;
        double cost;
        double partial_ObjFunc;
    };
    vector<move> moves;
    int feasibleMoves;

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

    void convertertoSolution(vector<double>& stat);

    void randomMutation();

    void getStatSolution(vector<double>& stat);

    void printProblem();

    void writeKPI(string path, string nameInstance, vector<double> stat);

    void writeSolution(string path);

    eFeasibleState isFeasible(string path);
};

#endif //COIOTE_HEURISTIC_HEURISTIC_H
