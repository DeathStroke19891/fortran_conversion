import numpy as np


def fermi(x):
    x =np.float128(1.0 / (np.exp(x) + 1.0))
    #print(x)
    return x


def sigdc_T(Z, W, tx, ty, tz, om_ref, Ef, Temp):
    N, nx, ny, nz = Z.shape[0], 22, 22, 5
    N_plane, N_line = nx * ny, nx

    Pi = np.arccos(-1.0)
    dsigz = np.zeros(10, dtype=np.complex128)
    om = np.zeros(10, dtype=np.float64)
    avsig = np.zeros(1000, dtype=np.float64)

    for kkk in range(1, 11):
        omega = om_ref * kkk
        dsigz[kkk - 1] = 0.0
        icounter = 0

        for mm in range(1, N):
            for nn in range(mm + 1, N + 1):
                sum1 = sum2 = sum3 = sum4 = sum5 = sum6 = 0.0

                if abs(W[mm - 1] - W[nn - 1]) >= 1e-9:
                    if abs(W[mm - 1] - W[nn - 1]) <= omega:
                        x = (W[mm - 1] - Ef) / Temp
                        y = (W[nn - 1] - Ef) / Temp
                        icounter += 1
                        f1 = fermi(x)
                        f2 = fermi(y)

                        for k in range(1, nz):
                            for j in range(1, ny):
                                for i in range(1, nx):
                                    sum2 = np.conj(
                                        Z[
                                            (k - 1) * N_plane
                                            + (j - 1) * N_line
                                            + i
                                            - 1,
                                            mm - 1,
                                        ]
                                    )
                                    sum3 = Z[
                                        k * N_plane + (j - 1) * N_line + i - 1, nn - 1
                                    ]
                                    sum4 = np.conj(
                                        Z[
                                            k * N_plane + (j - 1) * N_line + i - 1,
                                            mm - 1,
                                        ]
                                    )
                                    sum5 = Z[
                                        (k - 1) * N_plane + (j - 1) * N_line + i - 1,
                                        nn - 1,
                                    ]
                                    sum6 = (sum2 * sum3) - (sum4 * sum5)
                                    sum1 += sum6

                        dsigz[kkk - 1] += (
                            (np.abs(sum1) ** 2) * (f1 - f2) / (W[nn - 1] - W[mm - 1])
                        )

        dsigz[kkk - 1] = tz**2 * dsigz[kkk - 1] / N

    sigmax1 = dsigz[0]
    sigmax2 = dsigz[1]
    sigmax3 = dsigz[2]
    sigmax4 = dsigz[3]
    sigmax5 = dsigz[4]
    sigmax6 = dsigz[5]
    sigmax7 = dsigz[6]
    sigmax8 = dsigz[7]
    sigmax9 = dsigz[8]
    sigmax10 = dsigz[9]

    dcsigma1 = (dsigz[1] - dsigz[0]) / (om[1] - om[0])
    dcsigma2 = (dsigz[2] - dsigz[1]) / (om[2] - om[1])
    dcsigma3 = (dsigz[3] - dsigz[2]) / (om[3] - om[2])
    dcsigma4 = (dsigz[4] - dsigz[3]) / (om[4] - om[3])
    dcsigma5 = (dsigz[5] - dsigz[4]) / (om[5] - om[4])
    dcsigma6 = (dsigz[6] - dsigz[5]) / (om[6] - om[5])
    dcsigma7 = (dsigz[7] - dsigz[6]) / (om[7] - om[6])
    dcsigma8 = (dsigz[8] - dsigz[7]) / (om[8] - om[7])
    dcsigma9 = (dsigz[9] - dsigz[8]) / (om[9] - om[8])

    return (
        sigmax1,
        sigmax2,
        sigmax3,
        sigmax4,
        sigmax5,
        sigmax6,
        sigmax7,
        sigmax8,
        sigmax9,
        sigmax10,
        dcsigma1,
        dcsigma2,
        dcsigma3,
        dcsigma4,
        dcsigma5,
        dcsigma6,
        dcsigma7,
        dcsigma8,
        dcsigma9,
    )


# Define the input parameters
nx, ny, nz = 22, 22, 5
N = nx * ny * nz
Z = np.random.rand(N, N) + 1j * np.random.rand(N, N)
W = np.random.rand(N)
tx, ty, tz = 1.0, 1.0, 1.0
om_ref = 0.56
Ef = 0.0
Temp = 0.001

# Call the sigdc_T subroutine
results = sigdc_T(Z, W, tx, ty, tz, om_ref, Ef, Temp)

# Print the results
for i, result in enumerate(results):
    print(f"Result {i + 1}: {result}")
