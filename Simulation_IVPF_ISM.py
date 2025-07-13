import numpy as np
import matplotlib.pyplot as plt
import random

def safe_power(base, weights):
    base = np.maximum(base, 0)
    return np.prod(base ** weights)

def aggregate_positive_membership(mu_L, mu_U, weights):
    return [1 - safe_power(1 - mu_L, weights), 1 - safe_power(1 - mu_U, weights)]

def aggregate_neutral_membership(mu_L, mu_U, eta_L, eta_U, weights):
    return [
        safe_power(1 - mu_L, weights) - safe_power(1 - mu_L - eta_L, weights),
        safe_power(1 - mu_U, weights) - safe_power(1 - mu_U - eta_U, weights)
    ]

def aggregate_negative_membership(mu_L, mu_U, eta_L, eta_U, nu_L, nu_U, weights):
    return [
        safe_power(1 - mu_L - eta_L, weights) - safe_power(1 - mu_L - eta_L - nu_L, weights),
        safe_power(1 - mu_U - eta_U, weights) - safe_power(1 - mu_U - eta_U - nu_U, weights)
    ]

def ivpfowia_operator(ivpf_array, weights):
    mu_L = np.array([v[0][0] for v in ivpf_array])
    mu_U = np.array([v[0][1] for v in ivpf_array])
    eta_L = np.array([v[1][0] for v in ivpf_array])
    eta_U = np.array([v[1][1] for v in ivpf_array])
    nu_L = np.array([v[2][0] for v in ivpf_array])
    nu_U = np.array([v[2][1] for v in ivpf_array])

    μ = aggregate_positive_membership(mu_L, mu_U, weights)
    η = aggregate_neutral_membership(mu_L, mu_U, eta_L, eta_U, weights)
    ν = aggregate_negative_membership(mu_L, mu_U, eta_L, eta_U, nu_L, nu_U, weights)

    toplam_L = μ[0] + η[0] + ν[0]
    toplam_U = μ[1] + η[1] + ν[1]

    if toplam_L > 1:
        oran = 1 / toplam_L
        μ[0] *= oran
        η[0] *= oran
        ν[0] *= oran

    if toplam_U > 1:
        oran = 1 / toplam_U
        μ[1] *= oran
        η[1] *= oran
        ν[1] *= oran

    return [μ, η, ν]

def score_ivpf(ivpf):
    μ, η, ν = ivpf
    return ((μ[0] + μ[1]) / 2) - ((η[0] + η[1]) / 2) - ((ν[0] + ν[1]) / 2)

def run_ivpf_fism_simulation(replications=200):
    DSS = []

    for _ in range(replications):
        z = random.randint(5, 20)
        x = random.randint(10, 30)
        FCrispDec = np.zeros((x, x))
        IVPFCrispDec = np.zeros((x, x))
        weights = np.full(z, 1 / z)

        for i in range(x):
            for j in range(x):
                Fl = Fm = Fr = 0
                ivpf_array = []

                for _ in range(z):
                    val = random.choice([0, 1, 2, 3, 4]) if i != j else 4

                    fz = {
                        0: {"Le": 0, "Mi": 0, "Ri": 0.25},
                        1: {"Le": 0, "Mi": 0.25, "Ri": 0.5},
                        2: {"Le": 0.25, "Mi": 0.5, "Ri": 0.75},
                        3: {"Le": 0.5, "Mi": 0.75, "Ri": 1.0},
                        4: {"Le": 0.75, "Mi": 1.0, "Ri": 1.0}
                    }[val]
                    Fl += fz["Le"]
                    Fm += fz["Mi"]
                    Fr += fz["Ri"]

                    ivpf = {
                        0: [[0.01, 0.01], [0.40, 0.44], [0.50, 0.55]],
                        1: [[0.15, 0.20], [0.30, 0.35], [0.40, 0.45]],
                        2: [[0.35, 0.40], [0.20, 0.25], [0.30, 0.35]],
                        3: [[0.55, 0.60], [0.10, 0.15], [0.20, 0.25]],
                        4: [[0.75, 0.80], [0.01, 0.05], [0.10, 0.15]]
                    }[val]
                    ivpf_array.append(ivpf)

                FCrispDec[i, j] = (Fl + 2 * Fm + Fr) / (4 * z)
                IVPFCrispDec[i, j] = score_ivpf(ivpfowia_operator(ivpf_array, weights))

        fuzzy_thresh = np.mean(FCrispDec)
        ivpf_thresh = np.mean(IVPFCrispDec)

        FIRM = (FCrispDec >= fuzzy_thresh).astype(int)
        IVPFRM = (IVPFCrispDec >= ivpf_thresh).astype(int)
        np.fill_diagonal(FIRM, 1)
        np.fill_diagonal(IVPFRM, 1)

        JaccA = np.sum(FIRM)
        JaccB = np.sum(IVPFRM)
        Joint = np.sum(np.logical_and(FIRM, IVPFRM))
        dss = 2 * Joint / (JaccA + JaccB)
        DSS.append(dss)

    return DSS

# Simülasyonu çalıştır ve görselleştir
dss_values = run_ivpf_fism_simulation(1000)

plt.figure(figsize=(8, 6))
plt.boxplot(dss_values, patch_artist=True, showmeans=True, meanline=True)
plt.title("Dice-Sørensen Similarity: IVPF-ISM vs. FISM")
plt.ylabel("DSS")
plt.xticks([1], ["Simulation"])
plt.text(1.1, np.mean(dss_values), f"Mean: {np.mean(dss_values):.2f}", color="red")
plt.text(1.1, np.mean(dss_values) - np.std(dss_values), f"SD: {np.std(dss_values):.2f}", color="green")
plt.text(1.1, min(dss_values), f"Replications: {len(dss_values)}", color="blue")
plt.grid(True)
plt.tight_layout()
plt.savefig("IVPF_FISM_DSS_Boxplot.png", dpi=300)
plt.savefig("IVPF_FISM_DSS_Boxplot.pdf", bbox_inches="tight")
plt.show()
