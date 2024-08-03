void CWENOFV::RunRK3(double deltaT)
{
    // Copy current solution to temporary storage
    for (int ei = m_startElemX; ei != m_endElemX; ei++)
        for (int ej = m_startElemY; ej != m_endElemY; ej++)
            for (int r = 0; r != m_varNum; ++r)
                m_Un[ei][ej].vector[r] = m_Uh[ei][ej].vector[r];

    for (int step = 0; step != 3; step++)
    {
        // Assemble the right-hand side (RHS) for the equations
        assembleRHS();

        double a, b, c;
        // Set coefficients for the third-order Runge-Kutta scheme based on the step
        switch (step)
        {
        case 0:
            a = 1.0;
            b = 0.0;
            c = 1.0;
            break;
        case 1:
            a = 3.0 / 4.0;
            b = 1.0 / 4.0;
            c = 1.0 / 4.0;
            break;
        case 2:
            a = 1.0 / 3.0;
            b = 2.0 / 3.0;
            c = 2.0 / 3.0;
            break;
        default:
            std::cout << "Error: Invalid RK3 step" << std::endl;
            std::cin.get();
            exit(1);
            break;
        }

        // Update solution vector using the RK3 coefficients
        for (int ei = m_startElemX; ei != m_endElemX; ei++)
            for (int ej = m_startElemY; ej != m_endElemY; ej++)
                for (int r = 0; r < m_varNum; r++)
                    m_Uh[ei][ej].vector[r] = a * m_Un[ei][ej].vector[r] + b * m_Uh[ei][ej].vector[r] + c * deltaT * m_rhs[ei][ej].vector[r];
    }
}