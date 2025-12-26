#include "gauss_newton.hpp"



void GaussNewton::setX0(BodyVector x)
{
    x0 = x;
    x_best = x;
}


void GaussNewton::setJD(double jd)
{
    jd_start = jd;
}

BodyVector GaussNewton::getXbest() const
{
    return x_best;
}
 
BodyVector GaussNewton::scaledToPhysical(const BodyVector &x)
{
    BodyVector phys;
    phys.x.x = x.x.x * SCALE[0];
    phys.x.y = x.x.y * SCALE[1];
    phys.x.z = x.x.z * SCALE[2];

    phys.v.x = x.v.x * SCALE[3];
    phys.v.y = x.v.y * SCALE[4];
    phys.v.z = x.v.z * SCALE[5];
    return phys;
}

BodyVector GaussNewton::physicalToScaled(const BodyVector &x)
{
    BodyVector scal;
    scal.x.x = x.x.x / SCALE[0];
    scal.x.y = x.x.y / SCALE[1];
    scal.x.z = x.x.z / SCALE[2];

    scal.v.x = x.v.x / SCALE[3];
    scal.v.y = x.v.y / SCALE[4];
    scal.v.z = x.v.z / SCALE[5];
    return scal;
}



std::vector<double> GaussNewton::computeResiduals(const BodyVector &x)
{   
    BodyVector cur = scaledToPhysical(x);

    r.integrateOrbit(cur, jd_start, 2 * 365, 3600);

    std::vector<double> residuals;
    residuals.reserve(2 * a.size());

    for (auto [cur_jd_utc, Ra, Dec, code] : a)
    {
        double cur_jd_tdb = obs.convertUtcToTdb(cur_jd_utc, code);
        auto comp_pos_obs  = obs.getConvertObsForJD(cur_jd_tdb, code);
    
        auto Earth = eph.getStateVector((int)Ephemeris::CelestialBody::Earth, (int) Ephemeris::CelestialBody::SSB, cur_jd_tdb);

        auto Sun = eph.getStateVector((int)Ephemeris::CelestialBody::Sun, (int) Ephemeris::CelestialBody::SSB, cur_jd_tdb);

        auto ltc = r.applyAstrometricCorrections(cur_jd_tdb, comp_pos_obs, Earth, Sun);
        double Ra1, Dec1, dist1;
        r.cartToRaDec(ltc, Ra1, Dec1, dist1);
        residuals.push_back(Ra - Ra1);
        residuals.push_back(Dec - Dec1);
    }
    return residuals;
}

std::vector<std::vector<double>> GaussNewton::computeJ(const BodyVector &x)
{
    BodyVector cur = physicalToScaled(x);
    const double h = 1e-3;
    std::size_t n_res = 2 * a.size();
    std::size_t n_par = 6;

    std::vector<std::vector<double>> J(n_res, std::vector<double>(n_par, 0.0));
    
    for (std::size_t idx = 0; idx < n_par; ++idx)
    {
        BodyVector state_plus = cur;
        BodyVector state_minus = cur;
        
        switch(idx) 
        {
            case 0: state_plus.x.x += h; state_minus.x.x -= h; break;
            case 1: state_plus.x.y += h; state_minus.x.y -= h; break;
            case 2: state_plus.x.z += h; state_minus.x.z -= h; break;
            case 3: state_plus.v.x += h; state_minus.v.x -= h; break;
            case 4: state_plus.v.y += h; state_minus.v.y -= h; break;
            case 5: state_plus.v.z += h; state_minus.v.z -= h; break;
        }
        
        std::vector<double> r_plus = computeResiduals(state_plus);
        std::vector<double> r_minus = computeResiduals(state_minus);
        
        for (std::size_t i = 0; i < n_res; ++i) 
        {
            J[i][idx] = (r_plus[i] - r_minus[i]) / (2.0 * h);
        }    
    }
    
    return J;
}

bool GaussNewton::solveSystem(const std::vector<std::vector<double>> &A, 
                     const std::vector<double> &b, 
                     std::vector<double> &x)
{
    int n = A.size();
    x.resize(n, 0.0);
    
    // A = L * L^T
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));
    
    for (int i = 0; i < n; ++i) 
    {
        for (int j = 0; j <= i; ++j) 
        {
            double sum = 0.0;
            
            for (int k = 0; k < j; ++k) 
            {
                sum += L[i][k] * L[j][k];
            }
            
            if (i == j) {
                double diag = A[i][i] - sum;
                if (diag <= 0.0) 
                {
                    return false;
                }
                L[i][i] = sqrt(diag);
            } 
            else 
            {
                L[i][j] = (A[i][j] - sum) / L[j][j];
            }
        }
    }
    
    // L * z = b
    std::vector<double> z(n, 0.0);
    for (int i = 0; i < n; ++i) 
    {
        double sum = 0.0;
        for (int j = 0; j < i; ++j) 
        {
            sum += L[i][j] * z[j];
        }
        z[i] = (b[i] - sum) / L[i][i];
    }
    
    // L^T * x = z
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += L[j][i] * x[j];
        }
        x[i] = (z[i] - sum) / L[i][i];
    }
    
    return true;
}


int GaussNewton::start(std::size_t max_iters, double eps)
{
    BodyVector cur = x0;
    BodyVector best = cur;
    double best_rmse = 1e30;
    std::size_t n_obs = a.size();
    


    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Starting GaussNewton:\n";

    for (std::size_t iter = 0; iter < max_iters; ++iter)
    {
        std::cout << "═══════════════════════════════════════════════════════════\n";
        std::cout << "Iteration " << iter + 1 << "\n";
        std::vector<double> r = computeResiduals(physicalToScaled(cur));
        double mse = 0.0;
        for (auto now : r)
        {
            mse += now * now;
        }
        double rmse = sqrt(mse);
        std::cout << "RMSE " << rmse << "\n";
        if (sqrt(rmse) < best_rmse)
        {
            best_rmse = sqrt(rmse);
            best = cur;
        }
        auto J = computeJ(cur);
       
        int n_params = 6;
        std::vector<std::vector<double>> A(n_params, std::vector<double>(n_params, 0.0));
        std::vector<double> b(n_params, 0.0);
       
        // A = J^T * J
        for (int i = 0; i < n_params; ++i) 
        {
            for (int j = 0; j < n_params; ++j) 
            {
                double sum = 0.0;
                for (size_t k = 0; k < J.size(); ++k)
                {
                    sum += J[k][i] * J[k][j];
                }
                A[i][j] = sum;
            }
        }
       
        // b = J^T * r
        for (int i = 0; i < n_params; ++i) {
            double sum = 0.0;
            for (size_t k = 0; k < J.size(); ++k) {
                sum += J[k][i] * r[k];
            }
            b[i] = sum;
        }

        std::vector<double> delta_scaled(6, 0.0);
        bool solved = solveSystem(A, b, delta_scaled);
        
        if (!solved) 
        {
            std::cout << "⚠ Linear solver failed\n";
            return 0;
        }
        BodyVector delta_phys;
        delta_phys.x.x = delta_scaled[0] * SCALE[0];
        delta_phys.x.y = delta_scaled[1] * SCALE[1];
        delta_phys.x.z = delta_scaled[2] * SCALE[2];
        delta_phys.v.x = delta_scaled[3] * SCALE[3];
        delta_phys.v.y = delta_scaled[4] * SCALE[4];
        delta_phys.v.z = delta_scaled[5] * SCALE[5];

        cur.x.x -= delta_phys.x.x;
        cur.x.y -= delta_phys.x.y;
        cur.x.z -= delta_phys.x.z;
        cur.v.x -= delta_phys.v.x;
        cur.v.y -= delta_phys.v.y;
        cur.v.z -= delta_phys.v.z;

        double delta_pos_norm = sqrt(delta_phys.x.x * delta_phys.x.x +
                                     delta_phys.x.y * delta_phys.x.y +
                                     delta_phys.x.z * delta_phys.x.z);
        
        double delta_vel_norm = sqrt(delta_phys.v.x * delta_phys.v.x +
                                     delta_phys.v.y * delta_phys.v.y +
                                     delta_phys.v.z * delta_phys.v.z);
        
        double delta_total_norm = sqrt(delta_pos_norm * delta_pos_norm +
                                       delta_vel_norm * delta_vel_norm);
        
        std::cout << "Δ position norm: " << delta_pos_norm << " m\n";
        std::cout << "Δ velocity norm: " << delta_vel_norm << " m/s\n";
        std::cout << "Δ total norm:    " << delta_total_norm << "\n";
        std::cout << "Position: (" << cur.x.x << ", "
                  << cur.x.y << ", " << cur.x.z << ")\n";
        std::cout << "Velocity: (" << cur.v.x << ", "
                  << cur.v.y << ", " << cur.v.z << ")\n";
        if (delta_total_norm < 1e-3)
        {
            std::cout << "The vector changes by less than 1e-3 m\n";
            break;
        }
    }
    x_best = cur;
    return 1;
}

