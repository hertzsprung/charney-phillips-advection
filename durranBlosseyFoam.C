#include "fvCFD.H"
#include "simpleMatrix.H"

typedef struct
{
    const dimensionedScalar& x_min;
    const dimensionedScalar& x_max;
    const dimensionedScalar& z_min;
    const dimensionedScalar& z_max;
    const int nx;
    const int nz;
    const dimensionedScalar& dx;
    const dimensionedScalar& dz;
} geometry;

void initialise(Field< Field<scalar> >& phi, const geometry& mesh);
void write(Field< Field<scalar> >& phi, dimensionedScalar t, const geometry& mesh);

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #define dt runTime.deltaT()

    geometry mesh = {
       dimensionedScalar("x_min", dimLength, -150e3),
       dimensionedScalar("x_min", dimLength, 150e3),
       dimensionedScalar("z_min", dimLength, 0),
       dimensionedScalar("z_max", dimLength, 25e3),
       300,
       50,
       dimensionedScalar("dx", (mesh.x_max - mesh.x_min)/mesh.nx),
       dimensionedScalar("dz", (mesh.z_max - mesh.z_min)/mesh.nz),
    };

    const dimensionedScalar u0("u0", dimVelocity, 10);
    const dimensionedScalar cx("c", u0 * dt / mesh.dx);

    Info << "# t x theta" << endl;
   
    Field< Field<scalar> > theta_old(mesh.nz+1);
    Field< Field<scalar> > theta(mesh.nz+1);
    Field< Field<scalar> > theta_new(mesh.nx);

    dimensionedScalar t("t", dimTime, 0);
    initialise(theta, mesh);
    initialise(theta_old, mesh);
    initialise(theta_new, mesh);
    write(theta, t, mesh);

//    theta_old = theta;
//    theta_new = theta;
    scalar dt_multiplier = 0.5; // forward-in-time for the first timestep

    while (runTime.loop())
    {
        forAll(theta_new, K)
        {
            forAll(theta_new[K], I)
            {
                theta_new[K][I] = theta_old[K][I] - dt_multiplier*cx.value()*(theta[K][(I+1)%mesh.nx] - theta[0][(I-1)%mesh.nx]);
            }
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

void initialise(Field< Field<scalar> >& phi, const geometry& mesh)
{
    dimensionedScalar x0("x0", dimLength, -50e3);
    dimensionedScalar z0("z0", dimLength, 9e3);
    dimensionedScalar Ax("Ax", dimLength, 25e3);
    dimensionedScalar Az("Az", dimLength, 3e3);

    dimensionedScalar z = mesh.z_min;
    forAll(phi, K)
    {
        dimensionedScalar x = mesh.x_min + mesh.dx/2;
        phi[K] = Field<scalar>(mesh.nx); // how is this allocated?  will it always be in scope?  I need to learn me some more C++ :(
        forAll(phi[K], I)
        {
            dimensionedScalar r = sqrt(sqr((x-x0)/Ax) + sqr((z-z0)/Az));
            phi[K][I] = (r.value() <= 1) ? sqr(Foam::cos(M_PI*r.value()/2)) : 0;
            x += mesh.dx;
        }
        z += mesh.dz;
    }
}

void write(Field< Field<scalar> >& phi, dimensionedScalar t, const geometry& mesh)
{
    dimensionedScalar z = mesh.z_min;
    forAll(phi, K)
    {
        dimensionedScalar x = mesh.x_min + mesh.dx/2;
        forAll(phi[K], I)
        {
            Info << t.value() << " ";
            Info << x.value() << " " << z.value() << " " << phi[K][I] << endl;
            x += mesh.dx;
        }
        z += mesh.dz;
    }
    Info << endl << endl;
}
