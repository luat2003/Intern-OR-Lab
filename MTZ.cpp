// First formulation : MTZ formulation
// Miller–Tucker–Zemlin formulation
// BOOLEAN VARIABLE x[i][j]    For all i,j = 0, ..., n - 1
    //    x[i][j] == 1            If arc (i,j) is selected
    //    x[i][j] == 0            Otherwise
    //
// INTEGER VARIABLE u[i]      For all i = 0, ..., n - 1
    //    u[i] == k               Iff node i is the k-th node in the tour
    //    u[0] == 1
    //    u[i] in [2, ..., n]     For all i = 1, ... n - 1
    //
// OBJECTIVE FUNCTION
    //    Minimize sum((i,j), c[i][j] * x[i][j])
    //
// CONSTRAINTS
    //    1) sum(j, x[j][i]) == 1                    For all i
    //    2) sum(j, x[i][j]) == 1                    For all i
    //    3) u[i] - u[j] + n * x[i][j] <= n - 1   For all i,j = 1, ..., n - 1
    //       



#include<ilcplex\ilocplex.h>
#include<iostream>
ILOSTLBEGIN

using namespace std;
typedef IloArray<IloNumVarArray> IloNumVarArray2;
float distance(float x, float y) {
    return sqrt(x * x + y * y);
}

int main() {
    // input
    freopen("Test.txt", "r", stdin);


    int t, n = 0;
    float* coor_x = new float[100005];
    float* coor_y = new float[100005];
    // type 1
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
        IloNumVarArray u(env, n, 1, n, ILOINT);
        for (int i = 1; i < n; i++) {
            for (int j = 1; j < n; j++) {
                con3.add(u[i] - u[j] + n * x[i][j] <= n - 1);
            }
            con3.add(u[i] >= 2);
        }
        model.add(con3);
        model.add(u[0] == 1);

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

        IloNumArray rank(env, n);
        cplex.getValues(rank, u);
        env.out() << "Rank: " << rank << "\n";
        env.out() << "Tour by Rank: " << 1;
        int* res = new int[n];
        for (int i = 1; i < n; i++) {
            res[(int)rank[i] - 1] = i + 1;
        }
        for (int i = 1; i < n; i++) {
            env.out() << " -> " << res[i];
        }
        env.out() << "\n";


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