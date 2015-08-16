#include "fvCFD.H"
#include "simpleMatrix.H"

int main(int argc, char *argv[])
{
    simpleMatrix<scalar> m(4, 0, 0);
    Field<scalar>& rhs = m.source();
    rhs[0] = 0;
    rhs[1] = 2;
    rhs[2] = 4;
    rhs[3] = 6;

    for (int i=0; i<4; i++) {
        m[i][i] = 2;
    }

    const Field<scalar>& solution = m.solve();

    Info << solution << endl;

    return 0;
}

