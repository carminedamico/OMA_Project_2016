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
#include <deque>

using namespace std;

struct Data {
    double**** costs;
    int* n;
    int* activities;
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
    int initialActivities;
    deque<Agent> usedAgents;
    Cell(int i, int n, double Obj) {
        this->i = i;
        this->activities = n;
        this->initialActivities = n;
        this->partialObjFunc = Obj;
    }
};

class Solution {
public:
    int cellsWithActivitiesLeft;
    vector<Cell> cells;
    double objFunc;
    int*** usersCell;
};

class Bee {
public:
	int iterationsWithoutImprovements = 0;
	int iterationsToRestore = 0;
	Solution solution;
	vector<Bee*> onlookers;
};

class Heuristic{
private:
    int nTimeSteps;
    int nCustomerTypes;
    int nCells;

    Data problem;

    bool hasSolution;
    double epsilon;
    Solution B;
	vector<Bee> bees;
    int**** solution;

    void copyDataStructure(Solution* Y, Solution* X);

public:
    /**
     * Default constructor
     * @return Heuristic object
     */
    Heuristic(){};
    double* bestSolution;
    bool bestSolutionKnown;

    /**
     * Constructor from external file
     * @param path path of the external file cotaining the instance of the problem
     * @return
     */
    Heuristic(string path);

    // stat contiene tempo in posizione 0 e objfunction in posizione 1
    void Metaheuristic(vector<double>& stat);

    void solveGreedy(Solution *R);

    void gentlemanGreedy(Solution *R);

    void gentlemanAgreement(Solution *R, Solution *S);

    void gentlemenClub(Solution *R, Solution *S, int n);

    void getStatSolution(vector<double>& stat);

    void printProblem();

    void writeKPI(string path, string nameInstance, vector<double> stat);

    void writeSolution(string path);

    eFeasibleState isFeasible(string path);
};

#endif //COIOTE_HEURISTIC_HEURISTIC_H
