# Orthogonal Time Frequency Space (OTFS) Wireless Transceiver Simulation

An end-to-end MATLAB implementation evaluating the performance of **Orthogonal Time Frequency Space (OTFS)** modulation in high-mobility wireless environments. This project benchmarks OTFS performance across realistic multipath fading channels (3GPP EVA and Nakagami-$m$) using a two-dimensional Minimum Mean Square Error (2D MMSE) equalizer.

---

## 📌 Project Overview

In high-speed vehicular ($5G/6G$) and low-Earth orbit (LEO) satellite communications, conventional OFDM waveforms suffer from severe Inter-Carrier Interference (ICI) due to high Doppler shifts. OTFS addresses this by mapping information symbols onto the **Delay-Doppler (DD) domain**, transforming a time-varying multipath channel into a localized, non-fading delay-doppler channel matrix.

This repository contains full end-to-end Monte Carlo simulation code that:
1. Simulates multi-stage signal generation, channel filtering, and 2D equalization.
2. Evaluates Bit Error Rate (BER) performance under severe fading conditions across varying SNR levels.

---

## 🛠️ Key Technical Features

* **Grid & Modulation Setup:** Configured with a $16 \times 16$ Delay-Doppler grid using 4-QAM mapping.
* **Symplectic ISFFT/SFFT Processing:** Converts symbols between Delay-Doppler and Time-Frequency domains using 2D Fourier transform chains.
* **Realistic Channel Modeling:**
  * **3GPP Extended Vehicular A (EVA):** 9-path power-delay profile simulating mobile vehicular propagation.
  * **Nakagami-$m$ Fading:** Parametric fading envelope modeling ($m=3.0$) to benchmark performance against baseline Rayleigh ($m=1$) fading.
* **2D MMSE Equalization:** Employs two-dimensional matrix equalization (`fft2` / `ifft2`) in the Delay-Doppler grid to simultaneously suppress Inter-Symbol Interference (ISI) and Inter-Carrier Interference (ICI).
* **Deterministic Monte Carlo Analysis:** Uses multi-frame averaging with controlled statistical variance to produce smooth logarithmic BER waterfall curves.

---

## 📂 Repository Structure

* **`OTFS_System_Analysis.m`**  
  Simulates a single-frame 12-stage pipeline tracking grid matrices, channel impulse responses, noise additions, and constellation recovery at every processing point.
* **`OTFS_BER_Comparison.m`**  
  Runs multi-frame Monte Carlo simulations over a wide SNR range to generate comparative BER vs. SNR logarithmic curves for 3GPP EVA and Nakagami-$m$ channels.

---

## 🚀 How to Run

1. Clone or download this repository to your local machine:
   ```bash
   git clone [https://github.com/Sabaree2208/OTFS-Wireless-Simulation.git](https://github.com/Sabaree2208/OTFS-Wireless-Simulation.git)
   
