#include "matrix.h"

vector<vector<double>> minor(vector<vector<double>> &mat, int param_i, int param_j)
{
    int n = mat.size(), m = mat[0].size();
    if(n == 1) return {{0}};
    vector<vector<double>> ans(n-1, vector<double>(m-1));
        for(int i = 0; i < n; i++)
        {
            if(i == param_i) continue;
            for(int j = 0; j < m; j++)
            {
                if(j == param_j) continue;
                ans[i - (i >= param_i)][j - (j >= param_j)] = mat[i][j];
            }
        }
    return ans;
}

double det(vector<vector<double>> &mat)
{
    int n = mat.size(), m = mat[0].size();

    if(n != m) return 0;
    if(n == 1) return mat[0][0];

    double ans = 0;

    for(int t = 0; t < m; t++)
    {
        auto mi = minor(mat, 0, t);
        ans += ((t%2 == 0) ? 1 : -1)*mat[0][t]*det(mi);
    }
    return ans;
}

vector<vector<double>> transpose(vector<vector<double>> &mat)
{
    int n = mat.size(), m = mat[0].size();
    vector<vector<double>> ans(m, vector<double>(n));
    for(int i = 0; i < n; i++)
        for(int j = 0; j < m; j++)
            ans[i][j] = mat[j][i];
    return ans;
}

vector<vector<double>> inverse(vector<vector<double>> &mat)
{
    int n = mat.size(), m = mat[0].size();
    double d = det(mat);
    if(!d) return {{}};

    vector<vector<double>> adj(n, vector<double>(m));
    // calculate adjoint
    for(int i = 0; i < n; i++)
        for(int j = 0; j < m; j++)
        {
            auto mi = minor(mat, i, j);
            adj[i][j] = (((i+j)%2 == 0) ? 1 : -1)*det(mi);
        }
    // transpose
    for(int i = 0; i < n; i++)
        for(int j = i + 1; j < m; j++)
            std::swap(adj[i][j], adj[j][i]);

    for(int i = 0; i < n; i++)
        for(int j = 0; j < m; j++)
            adj[i][j] /= d;
    return adj;
}

vector<double> solveSystem(vector<vector<double>> &mat, vector<double> &b)
{
    if(mat[0].size() != b.size()) return {};
    auto inv = inverse(mat);
    int n = inv.size(), m = inv[0].size();
    vector<double> ans(m);
    for(int i = 0; i < n; i++)
    {
        double x = 0;
        for(int j = 0; j < m; j++)
            x += inv[i][j]*b[j];
        ans[i] = x;
    }
    return ans;
}