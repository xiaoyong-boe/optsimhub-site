/**
 * Tolerance Analysis for Thin-Film Multilayer Coatings (TMM only)
 * Ported from Python: tolerance_analysis.py + transfer_matrix.py
 *
 * Provides:
 *   TMMCalc.tmm(pol, n_list, d_list, theta_deg, wl_nm) → {R, T, A}
 *   TMMCalc.sweep(n_list, d_list, theta_deg, wl_start, wl_end, n_pts) → {wls, R, T}
 *
 *   NormalSampler / UniformSampler — per-layer perturbation generators
 *
 *   TolSimulator(layers, n_inc, n_sub, theta_inc, ...)
 *     - runSingle(wl_nm, n_samples, seed) → {R_samples, R_mean, R_std, R_p5, R_p95, ...}
 *     - runEnvelope(wl_start, wl_end, n_pts, n_samples, seed) → {wls, R_nom, R_p5, R_p95}
 *
 * Layer format (same as UI):
 *   {Name, Type:"iso"|"uniaxial"|"biaxial", n_o, n_e, n_z, theta_fast, theta_xz, d}
 */

// ============================================================
// TMM Transfer Matrix Calculator
// ============================================================
const TMMCalc = {
    /**
     * Compute R, T, A for a given polarization and layer stack.
     *
     * @param {string} pol - "s" or "p"
     * @param {number[]} n_list - [n_inc, n_layer1, ..., n_layerN, n_sub]
     * @param {number[]} d_list - [0.0, d_layer1, ..., d_layerN, 0.0]
     * @param {number} theta_deg - incidence angle in degrees
     * @param {number} wl_nm - wavelength in nm
     * @returns {{R: number, T: number, A: number}}
     */
    tmm(pol, n_list, d_list, theta_deg, wl_nm) {
        const theta_rad = theta_deg * Math.PI / 180;
        const n_inc = n_list[0];
        const n_sub = n_list[n_list.length - 1];

        // Snell's law angles through stack
        const sin_theta_inc = (typeof n_inc === 'number') ? Math.sin(theta_rad) : 0;
        const cos_theta = [];
        for (let i = 0; i < n_list.length; i++) {
            const ni = (typeof n_list[i] === 'number') ? n_list[i] : n_list[i].real;
            let sin_t = sin_theta_inc / Math.max(ni, 1e-9);
            sin_t = Math.max(-1, Math.min(1, sin_t));
            cos_theta.push(Math.sqrt(Math.max(0, 1 - sin_t * sin_t)));
        }

        // Admittance for each layer
        // s-pol: η = n * cosθ,  p-pol: η = n / cosθ
        const eta = [];
        for (let i = 0; i < n_list.length; i++) {
            const ni = n_list[i];
            if (pol === 's') {
                eta.push(ni * cos_theta[i]);
            } else {
                eta.push(ni / cos_theta[i]);
            }
        }

        // Characteristic matrix method — M is stored as 2×4:
        //   Row 0: [M00_re, M00_im, M01_re, M01_im]
        //   Row 1: [M10_re, M10_im, M11_re, M11_im]
        let M = [[1, 0, 0, 0], [0, 0, 1, 0]];

        for (let i = 1; i < n_list.length - 1; i++) {
            const ni = n_list[i];
            const di = d_list[i];
            const ct = cos_theta[i];

            const delta = 2 * Math.PI * ni * di * ct / wl_nm;
            const cos_d = Math.cos(delta);
            const sin_d = Math.sin(delta);
            const eta_i = eta[i];

            // Layer matrix: [[cosδ, -i*sinδ/η], [-i*η*sinδ, cosδ]]
            // Stored as 2×4: Row0=[re(M00),im(M00),re(M01),im(M01)], Row1=[re(M10),im(M10),re(M11),im(M11)]
            const L = [[cos_d, 0, 0, -sin_d / eta_i],
                       [0, -eta_i * sin_d, cos_d, 0]];

            // Post-multiply: M = M_prev @ L  (right multiply, matches TMM fix)
            if (i === 1) {
                M = L;
            } else {
                // M @ L: complex 2×2 matrix multiplication
                const l = L, m = M;
                // C00 = M00*L00 + M01*L10
                const c00_re = m[0][0]*l[0][0] - m[0][1]*l[0][1] + m[0][2]*l[1][0] - m[0][3]*l[1][1];
                const c00_im = m[0][0]*l[0][1] + m[0][1]*l[0][0] + m[0][2]*l[1][1] + m[0][3]*l[1][0];
                // C01 = M00*L01 + M01*L11
                const c01_re = m[0][0]*l[0][2] - m[0][1]*l[0][3] + m[0][2]*l[1][2] - m[0][3]*l[1][3];
                const c01_im = m[0][0]*l[0][3] + m[0][1]*l[0][2] + m[0][2]*l[1][3] + m[0][3]*l[1][2];
                // C10 = M10*L00 + M11*L10
                const c10_re = m[1][0]*l[0][0] - m[1][1]*l[0][1] + m[1][2]*l[1][0] - m[1][3]*l[1][1];
                const c10_im = m[1][0]*l[0][1] + m[1][1]*l[0][0] + m[1][2]*l[1][1] + m[1][3]*l[1][0];
                // C11 = M10*L01 + M11*L11
                const c11_re = m[1][0]*l[0][2] - m[1][1]*l[0][3] + m[1][2]*l[1][2] - m[1][3]*l[1][3];
                const c11_im = m[1][0]*l[0][3] + m[1][1]*l[0][2] + m[1][2]*l[1][3] + m[1][3]*l[1][2];

                M = [[c00_re, c00_im, c01_re, c01_im],
                     [c10_re, c10_im, c11_re, c11_im]];
            }
        }

        // After all layers: B = M00 + η_sub * M01, C = M10 + η_sub * M11
        const eta_sub = eta[eta.length - 1];
        const B_re = M[0][0] + eta_sub * M[0][2];
        const B_im = M[0][1] + eta_sub * M[0][3];
        const C_re = M[1][0] + eta_sub * M[1][2];
        const C_im = M[1][1] + eta_sub * M[1][3];

        // r = (η_inc * B - C) / (η_inc * B + C)
        const eta_inc = eta[0];
        const num_re = eta_inc * B_re - C_re;
        const num_im = eta_inc * B_im - C_im;
        const den_re = eta_inc * B_re + C_re;
        const den_im = eta_inc * B_im + C_im;

        const r_num2 = num_re * num_re + num_im * num_im;
        const r_den2 = den_re * den_re + den_im * den_im;
        const R = r_num2 / Math.max(r_den2, 1e-20);

        // t = 2*η_inc / den
        const t_num2 = 4 * eta_inc * eta_inc;
        const T = (eta_sub / eta_inc) * t_num2 / Math.max(r_den2, 1e-20);

        const A = Math.max(0, 1 - R - T);

        return { R: Math.min(R, 1), T: Math.min(T, 1), A };
    },

    /**
     * Sweep wavelength range, return R and T arrays.
     * @returns {{wls: number[], R: number[], T: number[]}}
     */
    sweep(n_list, d_list, theta_deg, wl_start, wl_end, n_pts) {
        const wls = [];
        const R = [];
        const T = [];
        const step = (wl_end - wl_start) / (n_pts - 1);
        for (let i = 0; i < n_pts; i++) {
            const wl = wl_start + i * step;
            const rs = TMMCalc.tmm("s", n_list, d_list, theta_deg, wl);
            const rp = TMMCalc.tmm("p", n_list, d_list, theta_deg, wl);
            wls.push(wl);
            R.push(0.5 * (rs.R + rp.R));
            T.push(0.5 * (rs.T + rp.T));
        }
        return { wls, R, T };
    },

    /**
     * Single-wavelength R, T (unpolarized average).
     */
    rt(wl_nm, n_list, d_list, theta_deg) {
        const rs = TMMCalc.tmm("s", n_list, d_list, theta_deg, wl_nm);
        const rp = TMMCalc.tmm("p", n_list, d_list, theta_deg, wl_nm);
        return { R: 0.5 * (rs.R + rp.R), T: 0.5 * (rs.T + rp.T) };
    }
};

// ============================================================
// Samplers
// ============================================================

/** Normal/Gaussian sampler (Box-Muller) */
class NormalSampler {
    sample(mean, sigma, size, rng) {
        if (size !== undefined) {
            const result = new Float64Array(size);
            for (let i = 0; i < size; i++) {
                let u1, u2;
                do { u1 = rng(); } while (u1 <= 0);
                do { u2 = rng(); } while (u2 <= 0);
                const z = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
                result[i] = mean + z * sigma;
            }
            return result;
        }
        let u1, u2;
        do { u1 = rng(); } while (u1 <= 0);
        do { u2 = rng(); } while (u2 <= 0);
        return mean + Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2) * sigma;
    }
}

/** Uniform sampler */
class UniformSampler {
    sample(mean, halfRange, size, rng) {
        if (size !== undefined) {
            const result = new Float64Array(size);
            for (let i = 0; i < size; i++) {
                result[i] = mean + (2 * rng() - 1) * halfRange;
            }
            return result;
        }
        return mean + (2 * rng() - 1) * halfRange;
    }
}

// ============================================================
// Tolerance Simulator
// ============================================================
class TolSimulator {
    /**
     * @param {object[]} layers - UI layer format
     * @param {number} n_inc - incident medium index
     * @param {number} n_sub - substrate index
     * @param {number} theta_inc - incidence angle (degrees)
     * @param {number|number[]} d_tol_pct - thickness tolerance (%)
     * @param {number|number[]} n_tol_pct - n_o tolerance (%)
     * @param {number|number[]} n_e_tol_pct - n_e tolerance (%) for uniaxial/biaxial
     * @param {number|number[]} n_z_tol_pct - n_z tolerance (%) for biaxial (TMM mode unused, reserved for Berreman)
     */
    constructor(layers, n_inc, n_sub, theta_inc,
                d_tol_pct, n_tol_pct, n_e_tol_pct, n_z_tol_pct) {
        this.layers = layers;
        this.n_inc = n_inc;
        this.n_sub = n_sub;
        this.theta_inc = theta_inc;
        this.n_layers = layers.length;

        this.d_tol = this._resolve(d_tol_pct, 5.0);
        this.n_tol = this._resolve(n_tol_pct, 1.0);
        this.n_e_tol = n_e_tol_pct !== undefined ? this._resolve(n_e_tol_pct, 1.0) : [...this.n_tol];
        this.n_z_tol = n_z_tol_pct !== undefined ? this._resolve(n_z_tol_pct, 1.0) : [...this.n_tol];

        // Pre-compute nominal values
        this.d_nom = layers.map(l => l.d);
        this.n_o_nom = layers.map(l => l.n_o);
        this.n_e_nom = layers.map(l => l.n_e || l.n_o);
        this.is_biaxial = layers.map(l => (l.Type || 'iso') === 'biaxial');

        // Build TMM arrays
        const { nList, dList } = this._toTmmFormat();
        this.baseNList = nList;
        this.baseDList = dList;

        this.d_sampler = new NormalSampler();
        this.n_sampler = new NormalSampler();
    }

    _resolve(val, def) {
        if (Array.isArray(val)) return val;
        if (typeof val === 'number') return new Array(this.n_layers).fill(val);
        return new Array(this.n_layers).fill(def);
    }

    _toTmmFormat() {
        const nList = [this.n_inc];
        const dList = [0];
        for (const l of this.layers) {
            const n_o = l.n_o || 1.5;
            const n_e = l.n_e || n_o;
            const isBi = (l.Type || 'iso') === 'biaxial';
            const n = (isBi && Math.abs(n_e - n_o) > 1e-6) ? Math.sqrt(n_o * n_e) : n_o;
            nList.push(n);
            dList.push(l.d);
        }
        nList.push(this.n_sub);
        dList.push(0);
        return { nList, dList };
    }

    /**
     * Run Monte Carlo at a single wavelength.
     * @returns object with R_samples, R_mean, R_std, R_p5, R_p95, R_min, R_max, ...
     */
    runSingle(wl_nm, n_samples, seed) {
        let rng = this._makeRng(seed);
        const R_samples = new Float64Array(n_samples);

        for (let s = 0; s < n_samples; s++) {
            const nList = [...this.baseNList];
            const dList = [...this.baseDList];

            for (let i = 0; i < this.n_layers; i++) {
                // Perturb thickness
                const d_scale = this.d_nom[i] * this.d_tol[i] / 300;
                const d_new = this.d_sampler.sample(this.d_nom[i], Math.max(d_scale, 1e-9), undefined, rng);
                dList[i + 1] = Math.max(1, d_new);

                // Perturb n
                const n_scale = this.n_o_nom[i] * this.n_tol[i] / 300;
                const no_new = this.n_sampler.sample(this.n_o_nom[i], Math.max(n_scale, 1e-9), undefined, rng);

                if (this.is_biaxial[i]) {
                    const ne_scale = this.n_e_nom[i] * this.n_e_tol[i] / 300;
                    const ne_new = this.n_sampler.sample(this.n_e_nom[i], Math.max(ne_scale, 1e-9), undefined, rng);
                    nList[i + 1] = Math.max(0.01, Math.sqrt(no_new * ne_new));
                } else {
                    nList[i + 1] = Math.max(0.01, no_new);
                }
            }

            const rt = TMMCalc.rt(wl_nm, nList, dList, this.theta_inc);
            R_samples[s] = rt.R;
        }

        const sorted = Array.from(R_samples).sort((a, b) => a - b);
        const idx5 = Math.floor(n_samples * 0.05);
        const idx95 = Math.floor(n_samples * 0.95);

        let sum = 0, sumSq = 0;
        for (let s = 0; s < n_samples; s++) {
            sum += R_samples[s];
            sumSq += R_samples[s] * R_samples[s];
        }
        const mean = sum / n_samples;
        const std = Math.sqrt(Math.max(0, sumSq / n_samples - mean * mean));

        // Yield: fraction of samples within ±10% of nominal R
        const R0 = TMMCalc.rt(wl_nm, [...this.baseNList], [...this.baseDList], this.theta_inc).R;
        const tolMargin = Math.max(0.01, R0 * 0.10); // at least 1% absolute margin
        const yieldCnt = R_samples.filter(r => r >= R0 - tolMargin && r <= R0 + tolMargin).length;

        return {
            R_samples, R_mean: mean, R_std: std,
            R_p5: sorted[idx5], R_p95: sorted[idx95],
            R_min: sorted[0], R_max: sorted[n_samples - 1],
            yield99: yieldCnt / n_samples,
            yield95: 0,
        };
    }

    /**
     * Run tolerance envelope over wavelength range.
     * @returns { wls, R_nom, R_p5, R_p95 }
     */
    runEnvelope(wl_start, wl_end, n_pts, n_samples, seed) {
        const rng = this._makeRng(seed || 42);
        const wls = [];
        const step = (wl_end - wl_start) / (n_pts - 1);
        for (let i = 0; i < n_pts; i++) wls.push(wl_start + i * step);

        const R_nom = new Float64Array(n_pts);
        const R_p5 = new Float64Array(n_pts);
        const R_p95 = new Float64Array(n_pts);

        for (let j = 0; j < n_pts; j++) {
            // Nominal
            const nomRt = TMMCalc.rt(wls[j], this.baseNList, this.baseDList, this.theta_inc);
            R_nom[j] = nomRt.R;

            // MC samples
            const vals = new Float64Array(n_samples);
            for (let s = 0; s < n_samples; s++) {
                const nList = [...this.baseNList];
                const dList = [...this.baseDList];
                for (let i = 0; i < this.n_layers; i++) {
                    const d_scale = this.d_nom[i] * this.d_tol[i] / 300;
                    dList[i + 1] = Math.max(1,
                        this.d_sampler.sample(this.d_nom[i], Math.max(d_scale, 1e-9), undefined, rng));
                    const n_scale = this.n_o_nom[i] * this.n_tol[i] / 300;
                    const no_new = this.n_sampler.sample(this.n_o_nom[i], Math.max(n_scale, 1e-9), undefined, rng);
                    if (this.is_biaxial[i]) {
                        const ne_scale = this.n_e_nom[i] * this.n_e_tol[i] / 300;
                        const ne_new = this.n_sampler.sample(this.n_e_nom[i], Math.max(ne_scale, 1e-9), undefined, rng);
                        nList[i + 1] = Math.max(0.01, Math.sqrt(no_new * ne_new));
                    } else {
                        nList[i + 1] = Math.max(0.01, no_new);
                    }
                }
                vals[s] = TMMCalc.rt(wls[j], nList, dList, this.theta_inc).R;
            }
            const sorted = Array.from(vals).sort((a, b) => a - b);
            R_p5[j] = sorted[Math.floor(n_samples * 0.05)];
            R_p95[j] = sorted[Math.floor(n_samples * 0.95)];
        }

        return { wls, R_nom, R_p5, R_p95 };
    }

    _makeRng(seed) {
        // Simple xoshiro128** PRNG
        let s0 = (seed || 42) >>> 0;
        let s1 = (s0 * 1812433253 + 1) >>> 0;
        let s2 = (s1 * 1812433253 + 1) >>> 0;
        let s3 = (s2 * 1812433253 + 1) >>> 0;
        return () => {
            const result = (Math.imul(s1 * 5, 1) << 7 | (s1 * 5) >>> 25) * 9;
            const t = s1 << 9;
            s2 ^= s0;
            s3 ^= s1;
            s1 ^= s2;
            s0 ^= s3;
            s2 ^= t;
            s3 = (s3 << 11 | s3 >>> 21);
            return (result >>> 0) / 4294967296;
        };
    }
}
