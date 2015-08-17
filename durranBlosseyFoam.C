#include "fvCFD.H"
#include "simpleMatrix.H"

typedef struct
{
    const dimensionedScalar& x_min;
    const dimensionedScalar& x_max;
    const int nx;
    const dimensionedScalar& dx;
} geometry;

void initialise(Field<scalar>& phi, const geometry& mesh);
void write(Field<scalar>& phi, dimensionedScalar t, const geometry& mesh);

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #define dt runTime.deltaT()

    geometry mesh = {
       dimensionedScalar("x_min", dimLength, -150e3),
       dimensionedScalar("x_min", dimLength, 150e3),
       300,
       dimensionedScalar("dx", (mesh.x_max - mesh.x_min)/mesh.nx)
    };
    
    Info << mesh.dx << endl;

    const dimensionedScalar u0("u0", dimVelocity, 10);
    const dimensionedScalar c("c", u0 * dt / mesh.dx);

    Info << "# Courant no " << c.value() << endl;
    Info << "# t x theta" << endl;
   
    Field<scalar> theta_old(mesh.nx);
    Field<scalar> theta(mesh.nx);
    Field<scalar> theta_new(mesh.nx);

    dimensionedScalar t("t", dimTime, 0);
    initialise(theta, mesh);
    write(theta, t, mesh);

    theta_old = theta;
    scalar dt_multiplier = 0.5; // forward-in-time for the first timestep

    while (runTime.loop())
    {
        forAll(theta_new, I)
        {
            theta_new[I] = theta_old[I] - dt_multiplier*c.value()*(theta[(I+1)%mesh.nx] - theta[(I-1)%mesh.nx]);
        }

        t += dt;
        dt_multiplier = 1;

        if (static_cast<int>(t.value()) % 5000 == 0)
        {
            write(theta_new, t, mesh);
        }

        theta_old = theta;
        theta = theta_new;
    }

    return 0;
}

void initialise(Field<scalar>& phi, const geometry& mesh)
{
    dimensionedScalar x0("x0", dimLength, -50e3);
    dimensionedScalar Ax("Ax", dimLength, 25e3);
//    scalar z0 = 9e3;
//    scalar Az = 3e3;

    dimensionedScalar x = mesh.x_min + mesh.dx/2;
    forAll(phi, I)
    {
        dimensionedScalar r = sqrt(sqr((x-x0)/Ax));
        phi[I] = (r.value() <= 1) ? sqr(Foam::cos(M_PI*r.value()/2)) : 0;
        x += mesh.dx;
    }
}

void write(Field<scalar>& phi, dimensionedScalar t, const geometry& mesh)
{
    dimensionedScalar x = mesh.x_min + mesh.dx/2;
    forAll(phi, I)
    {
        Info << t.value() << " ";
        Info << x.value() << " " << phi[I] << endl;
        x += mesh.dx;
    }
    Info << endl << endl;
}
