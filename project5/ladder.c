#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "matrix.h"

TMatrix * get_col_vec(TMatrix *x)
{
    int i = 0;
    int n = x->nrows;
    char buf[64];
    printf("\n\nPlease enter x0 as : x_0 x_1 x_2 etc.\n");
    fgets(buf, 64, stdin);
    char *token = strtok(buf, " ");
    do
    {
        // set x data
        x->data[i][0] = atoi(token);
        // get next token (x_n)
        token = strtok(NULL, " ");
        i++;
    } while(token != NULL && i <= n);
    // warning
    if(i != n)
    {
        printf("%d : WARNING the X vector has NOT been fully initialized, no garantee to the future behavior of the program...\n", i);
        exit(1);
    }
    return x;
}

void construct_a(TMatrix *a, double *g_s, double *g_g)
{
    // construct a
    int n = a->nrows;
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            if(i == 0)
            {
                if(j == 0)
                    a->data[i][j] = g_s[0];
                else if(j == 1)
                    a->data[i][j] = -g_s[0];
                else
                    a->data[i][j] = 0;
            }
            else if(i < n - 1)
            {
                if(j < i - 1)
                    a->data[i][j] = 0;
                else if(j == i - 1)
                    a->data[i][j] = g_s[i - 1];
                else if(j == i)
                    a->data[i][j] = -1* (g_s[i - 1] + g_s[i] + g_g[i - 1]);
                else if(j == i + 1)
                    a->data[i][j] = g_s[i];
            }
            else if(i == n - 1)
            {
                if(j < n - 2)
                    a->data[i][j] = 0;
                else if(j == n - 2)
                    a->data[i][j] = g_s[i - 1];
                else if(j ==  n - 1)
                    a->data[i][j] = -1 * (g_s[i - 1] + g_g[i - 1]);
            }
        }
    }
    // display a
    printf("\na is :\n\n");
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++) printf("%0.4f ", a->data[i][j]);
        printf("\n");
    }
}

double L2_norm_col_vec(TMatrix *x)
{
    double sum = 0;
    for(int i = 0; i < x->nrows; sum += pow(x->data[i][0], 2), i++);
    return sqrt(sum);
}

void mod_gauss_seidel(TMatrix *a, TMatrix *x, TMatrix *b, double tol)
{
    int n = a->nrows;
    TMatrix *L, *U, *l;
    L = newMatrix(n, n);
    U = newMatrix(n, n);
    // populate L and U
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            if(j > i)
            {
                L->data[i][j] = 0;
                U->data[i][j] = a->data[i][j];
            } else
            {
                L->data[i][j] = a->data[i][j];
                U->data[i][j] = 0;
            } 
        }
    }
    // display L
    printf("\nL is :\n\n");
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++) printf("%0.4f ", L->data[i][j]);
        printf("\n");
    }
    // find L inverse
    l = invertMatrix(L);
    // display L
    printf("\nl is :\n\n");
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++) printf("%0.4f ", l->data[i][j]);
        printf("\n");
    }
    // INVERSE FAILING, NEEDS WORK
}


TMatrix * gauss_seidel(TMatrix *a, TMatrix *x, TMatrix *b, double tol)
{
    int n = x->nrows;
    TMatrix *x_neg = newMatrix(n, 1);
    TMatrix *x_new = newMatrix(n, 1);
    TMatrix *dif = newMatrix(n, 1);
    // start x_new as x
    for(int i = 0; i < n; i++) x_new->data[i][0] = x->data[i][0];
    do
    {
        // update x
        for(int i = 0; i < n; i++) x->data[i][0] = x_new->data[i][0];
        // construct x_new
        for(int i = 0; i < n; i++)
        {
            double sum1, sum2 = 0;
            // forward substitution, assuming no zero diagonal entries.
            for(int j = 0; j < i - 1; j++) sum1 += a->data[i][j] * x_new->data[j][0];
            for(int j = i + 1; j < n; j++) sum2 += a->data[i][j] * x->data[j][0];
            x_new->data[i][0] = (b->data[i][0] - sum1 - sum2) / a->data[i][i];
        }
        // negate x, find interstep error
        for(int i = 0; i < n; i++) x_neg->data[i][0] = x->data[i][0] * -1;
        dif = addMatrix(x_new, x_neg);
    } while(L2_norm_col_vec(dif) > tol);
    // return last iteration
    return x_new;
}

int main(int argc, char **argv)
{
    // A in R nxn, Ax=b
    // get n
    int n;
    char buf[64];
    printf("Please enter n, the dimension of the space we are working in : ");
    fgets(buf, 64, stdin);
    n = atoi(buf);

    TMatrix *a, *b, *x, *y, *y_p;
    a = newMatrix(n, n);
    x = newMatrix(n, 1);
    b = newMatrix(n, 1);
    y = newMatrix(n, 1);
    int i = 0;

    // testing
    for(i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++) if(i == j) a->data[i][j] = i + 1; else a->data[i][j] = 0;
        x->data[i][0] = 0;
        b->data[i][0] = i + 1;
        y->data[i][0] = 1;
    }
    y = gauss_seidel(a, x, b, 0.1);
    printf("Should get 1's\n");
    for(i = 0; i < n; i++) printf("%0.4f ", y->data[i][0]);
    printf("\n");

    freeMatrix(a);
    freeMatrix(b);
    freeMatrix(x);
    freeMatrix(y);

    // ENDOF testing
    a = newMatrix(n, n);
    x = newMatrix(n, 1);
    b = newMatrix(n, 1);
    y = newMatrix(n, 1);
    i = 0;

    // circuit setup
    double r_s[n], r_g[n], g_s[n], g_g[n];
    r_s[0] = 400;
    r_g[0] = 50;
    g_s[0] = 1 / r_s[0];
    g_g[0] = 1/ r_g[0];
    for(i = 1; i < n; i++) 
    {
        r_s[i] = r_s[i-1] / 4;
        g_s[i] = 1 / r_s[i];
        r_g[i] = r_g[i-1] / 2;
        g_g[i] = 1 / r_g[i];
    }

    // CASE 1
    // get initial guess
    for(i = 0; i < n; i++) x->data[i][0] = 0;
    // construct b
    b->data[0][0] = 2;
    for(i = 1; i < n; i++) b->data[i][0] = 0;
    // display b
    printf("b is : ");
    for(i = 0; i < n; i++) printf("%0.4f ", b->data[i][0]);
    // construct a
    construct_a(a, g_s, g_g);
    // test gauss_seidel, negate answer to find difference
    y_p = gauss_seidel(a, x, b, 0.001);
    mod_gauss_seidel(a, x, b, 0.001);
    printf("\nx_n is :\n");
    for(i = 0; i < n; i++) printf("%0.4f ", y_p->data[i][0]);

    // CASE 2
    // get initial guess
    for(i = 0; i < n; i++) x->data[i][0] = 0;
    // construct b
    b->data[0][0] = 1;
    for(i = 1; i < n; i++) b->data[i][0] = 0;
    // display b
    printf("b is : ");
    for(i = 0; i < n; i++) printf("%0.4f ", b->data[i][0]);
    // reconstruct g_s, g_g
    for(i = 0; i < n; i++)
    {
        g_s[i] = 1;
        g_g[i] = 1;
    }
    // construct a
    construct_a(a, g_s, g_g);
    // test gauss_seidel, negate answer to find difference
    y_p = gauss_seidel(a, x, b, 0.001);
    mod_gauss_seidel(a, x, b, 0.001);
    printf("\nx_n is :\n");
    for(i = 0; i < n; i++) printf("%0.4f ", y_p->data[i][0]);

    // free memory
    freeMatrix(a);
    freeMatrix(b);
    freeMatrix(x);
    freeMatrix(y);
    freeMatrix(y_p);

    return 0;
}
