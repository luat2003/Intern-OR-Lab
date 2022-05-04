// Second formulation : DFJ formulation
// Dantzig–Fulkerson–Johnson formulation
// BOOLEAN VARIABLE x[i][j]    For all i,j = 0, ..., n - 1
    //    x[i][j] == 1            If arc (i,j) is selected
    //    x[i][j] == 0            Otherwise

// OBJECTIVE FUNCTION
    //    Minimize sum((i,j), c[i][j] * x[i][j])
    //
// CONSTRAINTS
    //    1) sum(j, x[j][i]) == 1                    For all i
    //    2) sum(j, x[i][j]) == 1                    For all i
    //    3) sum((i in Q,j in Q), x[i][j]) <= |Q| - 1 for all subset Q of {1,2,..,n} 
    //   1 < |Q| < n 


#include<ilcplex\ilocplex.h>
#include<iostream>
ILOSTLBEGIN

using namespace std;
typedef IloArray<IloNumVarArray> IloNumVarArray2;
float distance(float x, float y) {
    return sqrt(x * x + y * y);
}

int count_bitset(int n) {
    int cnt = 0;
    while (n) {
        if (n & 1) {
            cnt++;
        }
        n >>= 1;
    }
    return cnt;
}

int main() {
    // input
    freopen("Test.txt", "r", stdin);

    int t, n = 0;

    // type 1
    float* coor_x = new float[100005];
    float* coor_y = new float[100005];
    while (cin >> t) {
        cin >> coor_x[n] >> coor_y[n];
        n++;
    }
    float** c = new float* [n];
    for (int i = 0; i < n; i++) {
        c[i] = new float[n];
        for (int j = 0; j < n; j++) {
            c[i][j] = distance(coor_x[i] - coor_x[j], coor_y[i] - coor_y[j]);
        }
    }

    // type 2
    /*
    float** c = new float* [n];
    for (int i = 0; i < n; i++) {
        c[i] = new float[n];
        for (int j = 0; j < n; j++) {
            cin >> c[i][j];
        }
    }
    */

    IloEnv env;
    try {
        IloModel model(env);

        IloNumVarArray2 x(env, n);
        for (int i = 0; i < n; i++) {
            x[i] = IloNumVarArray(env, n, 0, 1, ILOBOOL);
        }

        IloExpr expr(env);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                expr += x[i][j] * c[i][j];
            }
        }
        model.add(IloMinimize(env, expr));
        expr.end();


        IloRangeArray con1(env);
        for (int i = 0; i < n; i++) {
            IloExpr exp0(env);
            for (int j = 0; j < n; j++) {
                exp0 += x[i][j];
            }
            con1.add(exp0 == 1);
        }
        model.add(con1);

        IloRangeArray con2(env);
        for (int j = 0; j < n; j++) {
            IloExpr exp0(env);
            for (int i = 0; i < n; i++) {
                exp0 += x[i][j];
            }
            con2.add(exp0 == 1);
        }
        model.add(con2);

        IloRangeArray con3(env);
        for (int i = 2; i < (1<<n)-1; i++) {
            int cnt = count_bitset(i);
            if (cnt > 1 && cnt < n) {
                IloExpr exp1(env);
                for (int l = 0; l < n; l++) {
                    if ( (1 << l) & i ) {
                        for (int j = 0; j < n; j++) {
                            if ( (1 << j) & i ) {
                                exp1 += x[l][j];
                            }
                        }
                    }
                }
                con3.add(exp1 <= cnt - 1);
            }
        }
        model.add(con3);

        IloCplex cplex(model); // create cplex object w/model
        cplex.solve(); // solve the model

        env.out() << "Solution status = " << cplex.getStatus() << "\n";
        env.out() << "Solution value  = " << cplex.getObjValue() << "\n";

        env.out() << "Matrix x = " << "\n";
        IloArray<IloNumArray> vals(env, n);
        for (int i = 0; i < n; i++) {
            vals[i] = IloNumArray(env, n, 0, 1, ILOBOOL);
            cplex.getValues(vals[i], x[i]);
            for (int j = 0; j < n; j++) {
                env.out() << abs(vals[i][j]) << " ";
            }
            env.out() << "\n";
        }

        IloInt now = 0, cnt = n;
        env.out() << "\n" << "Traveling by matrix x : " << 1;
        while (cnt--) {
            for (int i = 0; i < n; i++) {
                if (vals[now][i]) {
                    env.out() << " -> " << i + 1;
                    now = i;
                    break;
                }
            }
        }
        cout << "\n";




        cplex.exportModel("ipex1.lp");
    }
    catch (IloException& e) {
        cerr << "Concert exception caught: " << e << endl;
    }
    catch (...) {
        cerr << "Unknown exception caught" << endl;
    }

    env.end();

    return 0;
}