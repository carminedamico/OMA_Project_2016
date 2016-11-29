#include <list>
#include <iostream>

#include "utils.h"
#include "heuristic.h"

using namespace std;

int main(int argc,char *argv[]){
	bool _test = false;
	string _inPath;
	string _outPath;
	string _solPath;

	// Read input parameters
	for(int i=1; i< argc;i++) {
		if(strcmp(argv[i], "-i") == 0)
			_inPath = argv[++i];
		else if(strcmp(argv[i], "-o") == 0)
			_outPath = argv[++i];
		else if(strcmp(argv[i], "-s") == 0)
			_solPath = argv[++i];
		else if(strcmp(argv[i], "-test") == 0)
			_test = true;
	}

	if(_inPath.empty() || _outPath.empty() || _solPath.empty()) {
		cout << "------------------------------ " << endl;
		cout << "CMD" << endl;
		cout << "------------------------------ " << endl;
		cout << "-i path of the instance file or solution to test" << endl;
		cout << "-o path of the output file" << endl;
		cout << "-s path of the output solution file" << endl;
        return 1;
	}

	if(!_test) {
		// Read the instance file
		Heuristic _heuristic(_inPath);
		// Solve the problem
		vector<double> stat;
		_heuristic.Metaheuristic(stat);
        _heuristic.getStatSolution(stat);
		// Write KPI of solution
        string instanceName = splitpath(_inPath);
        _heuristic.writeKPI(_outPath, instanceName, stat);

        //
        cout << _inPath << ": \n";
        cout << "\tTIME: " << stat[1] << "\n";
        cout << "\tSOLUTION FOUNDED: " << stat[0] << "\n";
        cout << "\tBEST SOLUTION: " << _heuristic.bestSolution[0] << "\n";
        cout << "\tGAP: " << (double) ((stat[0] * 100) /  _heuristic.bestSolution[0]) - 100<< "%\n";
        for (int i = 2; i < stat.size(); i++) {
            int x = _heuristic.bestSolution[i - 1] - stat[i];
            if (x > 0) {
                cout << "\t+";
            } else cout << "\t";
            cout << x << " USERS OF TYPE " << i-1 << " USED\n";
        }

        _heuristic.writeSolution(_solPath);

        Heuristic _heuristicX = Heuristic(_inPath);
        // Read the solution file
        eFeasibleState _feasibility = _heuristicX.isFeasible(_solPath);
        switch(_feasibility) {
            case FEASIBLE:
                cout << "Solution is feasible" << endl;
                break;
            case NOT_FEASIBLE_DEMAND:
                cout << "Solution is not feasible: demand not satisfied" << endl;
                break;
            case NOT_FEASIBLE_FLOW:
                cout << "Solution is not feasible: flow constraint not satisfied" << endl;
                break;
            case NOT_FEASIBLE_USERS:
                cout << "Solution is not feasible: exceeded number of available users" << endl;
                break;
        }
	}
	return 0;
}
