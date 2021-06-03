#!/usr/bin/env python3
from functools import wraps
import matplotlib.pyplot as plt
from numpy import empty, exp, inf, isreal, linspace
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar
from time import perf_counter


class Katz(object):
    def __init__(self):
        self.A1 = 1085.7 + 273.15
        self.A2 = 132.9
        self.A3 = -5.1
        self.B1 = 1475 + 273.15
        self.B2 = 80
        self.B3 = -3.2
        self.C1 = 1780 + 273.15
        self.C2 = 45
        self.C3 = -2
        self.beta1 = 1.5
        self.beta2 = 1.5
        self.D_H2O = 0.01
        self.X_H2O_bulk = 0
        self.gam = 0.75
        self.K = 43
        self.ki1 = 12
        self.ki2 = 1
        self.lam = 0.6
        self.M_cpx = 0.17
        self.r0 = 0.5
        self.r1 = 0.08
        self.c_P = 1000  # Only for KatzPTF
        self.alpha_s = 4e-5  # Only for KatzPTF
        self.alpha_f = 6.8e-5  # Only for KatzPTF
        self.rho_s = 3300  # Only for KatzPTF
        self.rho_f = 2900  # Only for KatzPTF
        self.deltaS = 300  # Only for KatzPTF

    def updateConst(KatzFunc):
        @wraps(KatzFunc)
        def KatzFuncWrapper(*args, **kwargs):
            if kwargs.get('inputConst'):
                for key, value in kwargs['inputConst'].items():
                    args[0].__dict__[key] = value
                del kwargs['inputConst']
            return KatzFunc(*args, **kwargs)
        return KatzFuncWrapper

    def calcX_H2O(self, F):
        return self.X_H2O_bulk / (self.D_H2O + F * (1 - self.D_H2O))

    def calcSolLiqCpxOut(self, presGPa):
        T_sol = self.A1 + self.A2 * presGPa + self.A3 * presGPa ** 2
        T_liq_lherz = self.B1 + self.B2 * presGPa + self.B3 * presGPa ** 2
        T_liq = self.C1 + self.C2 * presGPa + self.C3 * presGPa ** 2
        F_cpx_out = self.M_cpx / (self.r0 + self.r1 * presGPa)
        T_cpx_out = (F_cpx_out ** (1 / self.beta1) * (T_liq_lherz - T_sol)
                     + T_sol)
        return T_sol, T_liq_lherz, T_liq, F_cpx_out, T_cpx_out

    def checkWaterSat(self, presGPa, temp, F, T_sol, T_liq_lherz):
        X_H2O_sat = self.ki1 * presGPa ** self.lam + self.ki2 * presGPa
        if self.calcX_H2O(F) > X_H2O_sat:
            tempPrime = ((temp - T_sol + self.K * X_H2O_sat ** self.gam)
                         / (T_liq_lherz - T_sol))
            F = tempPrime ** self.beta1 if tempPrime > 0 else 0
        return F

    @updateConst
    def KatzPT(self, presGPa, temp):
        def checkCpx(F):
            return ((temp - T_sol + self.K * self.calcX_H2O(F) ** self.gam)
                    / (T_liq_lherz - T_sol))

        def checkOpx(F):
            return ((temp - T_cpx_out + self.K * self.calcX_H2O(F) ** self.gam)
                    / (T_liq - T_cpx_out))

        def detBracket(start, stop, nbIncrmt, funcCheck, func):
            F = start
            for i in range(nbIncrmt):
                incrmt = 10 ** -(i + 1)
                funcCheckStart = funcCheck(start)
                if funcCheckStart <= 0:
                    return None
                funcStart = func(start, funcCheckStart)
                while F < stop:
                    funcCheckF = funcCheck(F)
                    if funcCheckF >= 0:
                        funcF = func(F, funcCheckF)
                        if isreal(funcF):
                            if funcF * funcStart < 0:
                                return [F - incrmt, F]
                            else:
                                F += incrmt
                        else:
                            break
                    else:
                        break
                F -= incrmt
            return [F, F + incrmt]

        def funcCpx(F, *args):
            checkCpxF = args[0] if args else checkCpx(F)
            return F - checkCpxF ** self.beta1

        def funcOpx(F, *args):
            checkOpxF = args[0] if args else checkOpx(F)
            return F - F_cpx_out - (1 - F_cpx_out) * checkOpxF ** self.beta2

        if presGPa > 8:
            return 0
        T_sol, T_liq_lherz, T_liq, F_cpx_out, T_cpx_out = (
            self.calcSolLiqCpxOut(presGPa))
        bracket = detBracket(0, F_cpx_out, 9, checkCpx, funcCpx)
        if bracket is None:
            return 0
        elif bracket[1] < F_cpx_out:
            F = root_scalar(funcCpx, method='brentq', bracket=bracket).root
            assert F >= 0
            if F > 0:
                F = self.checkWaterSat(presGPa, temp, F, T_sol, T_liq_lherz)
        else:
            bracket = detBracket(F_cpx_out, 1, 9, checkOpx, funcOpx)
            if bracket[1] >= 1:
                return 1
            F = root_scalar(funcOpx, method='brentq', bracket=bracket).root
        assert F >= 0 and F < 1
        return F

    @updateConst
    def KatzPTF(self, presGPaStart, presGPaEnd, tempStart, Fstart, dTdP_GPa):
        def deriv(t, y):
            presGPa = t
            temp, F = y
            T_sol, T_liq_lherz, T_liq, F_cpx_out, T_cpx_out = (
                self.calcSolLiqCpxOut(presGPa))
            if F <= 0:
                F = self.KatzPT(presGPa, temp)
            elif F < F_cpx_out:
                F = self.checkWaterSat(presGPa, temp, F, T_sol, T_liq_lherz)
            if F == 0:
                return dTdP_GPa, 0
            dT_sol = self.A2 + 2 * self.A3 * presGPa
            dT_liq_lherz = self.B2 + 2 * self.B3 * presGPa
            dT_liq = self.C2 + 2 * self.C3 * presGPa
            dH2O = (self.gam * self.K * self.X_H2O_bulk ** self.gam
                    * (1 - self.D_H2O)
                    / (self.D_H2O + F * (1 - self.D_H2O)) ** (self.gam + 1))
            if F < F_cpx_out:
                dTdP_F = (F ** (1 / self.beta1) * (dT_liq_lherz - dT_sol)
                          + dT_sol)
                dTdF_P = (F ** ((1 - self.beta1) / self.beta1)
                          * (T_liq_lherz - T_sol) / self.beta1 + dH2O)
            else:
                A = ((F - F_cpx_out) / (1 - F_cpx_out)) ** (1 / self.beta2)
                B = T_liq - T_cpx_out
                dAdP_F = (F_cpx_out ** 2 * self.r1 / self.M_cpx / self.beta2
                          * (F - F_cpx_out) ** (1 / self.beta2)
                          * (1 / (F - F_cpx_out) - 1 / (1 - F_cpx_out))
                          / (1 - F_cpx_out) ** (1 / self.beta2))
                dCdP_F = (-F_cpx_out ** ((self.beta1 + 1) / self.beta1)
                          * self.r1 / self.M_cpx / self.beta1
                          * (T_liq_lherz - T_sol)
                          + F_cpx_out ** (1 / self.beta1)
                          * (dT_liq_lherz - dT_sol) + dT_sol)
                dBdP_F = dT_liq - dCdP_F
                dTdP_F = dAdP_F * B + A * dBdP_F + dCdP_F
                dTdF_P = ((F - F_cpx_out) ** ((1 - self.beta2) / self.beta2)
                          / self.beta2 * (T_liq - T_cpx_out)
                          / (1 - F_cpx_out) ** (1 / self.beta2) + dH2O)
            dFdP_S = ((dTdP_GPa - dTdP_F)
                      / (temp * self.deltaS / self.c_P + dTdF_P))
            dTdP_S = dTdP_GPa - self.deltaS * dFdP_S * temp / self.c_P
            return dTdP_S, dFdP_S

        sol = solve_ivp(deriv, [presGPaStart, presGPaEnd], [tempStart, Fstart],
                        method='LSODA', first_step=None, max_step=inf,
                        t_eval=None, dense_output=True, atol=1e-6, rtol=1e-5)
        return sol.sol


class KatzFigures(object):
    @staticmethod
    def Figure2():
        temp = linspace(1000, 2000, 1001) + 273.15
        F = empty(temp.shape)
        katz = Katz()
        for pres in [0, 1, 2, 3]:
            begin = perf_counter()
            for i, T in enumerate(temp):
                F[i] = katz.KatzPT(pres, T)
            print(perf_counter() - begin)
            plt.plot(temp - 273.15, F, label=f'{pres} GPa')
        plt.legend(loc='upper left')
        plt.grid()
        plt.gca().set_xlim(temp[0] - 273.15, temp[-1] - 273.15)
        plt.gca().set_ylim(0, 1)
        plt.show()

    @staticmethod
    def Figure3():
        pres = linspace(8, 0, 401)
        temp = empty(pres.shape)
        katz = Katz()
        for water in [0, 0.05, 0.1, 0.3, 0.5]:
            begin = perf_counter()
            for i, P in enumerate(pres):
                T = 900 + 273.15
                while katz.KatzPT(P, T, inputConst={'X_H2O_bulk': water}) == 0:
                    T += 0.5
                temp[i] = T - 0.25
            print(perf_counter() - begin)
            plt.plot(temp - 273.15, pres, label=f'{water} bulk wt%')
        plt.legend(loc='upper right')
        plt.grid()
        plt.gca().set_xlim(1000 - 1.2 / 1.8 * 200, 1800 + 1.1 / 1.8 * 200)
        plt.gca().set_ylim(pres[0], pres[-1])
        plt.show()
        plt.show()

    @staticmethod
    def Figure4():
        temp = linspace(900, 1400, 501) + 273.15
        F = empty(temp.shape)
        katz = Katz()
        for water in [0, 0.02, 0.05, 0.1, 0.3]:
            begin = perf_counter()
            for i, T in enumerate(temp):
                F[i] = katz.KatzPT(1, T, inputConst={'X_H2O_bulk': water})
            print(perf_counter() - begin)
            plt.plot(temp - 273.15, F, label=f'{water} bulk wt%')
        plt.legend(loc='upper left')
        plt.grid()
        plt.gca().set_xlim(temp[0] - 273.15, temp[-1] - 273.15)
        plt.gca().set_ylim(0, 0.4)
        plt.show()

    @staticmethod
    def Figure11():
        fig, (ax, bx) = plt.subplots(nrows=1, ncols=2, sharey=True)
        presGPa = linspace(0, 6, 1001)
        katz = Katz()
        for poTemp, colour in zip([1250, 1350, 1450],
                                  ['tab:blue', 'tab:green', 'tab:orange']):
            temp = (poTemp + 273.15) * exp(katz.alpha_s * 6e9
                                           / katz.c_P / katz.rho_s)
            for water in linspace(0, 0.02, 4):
                begin = perf_counter()
                sol = katz.KatzPTF(6, 0, temp, 0, katz.alpha_s * temp
                                   / katz.rho_s / katz.c_P * 1e9,
                                   inputConst={'X_H2O_bulk': water,
                                               'M_cpx': 0.1})
                print(perf_counter() - begin)
                ax.plot(sol(presGPa)[1, :], presGPa, color=colour)
                bx.plot(sol(presGPa)[0, :] - 273.15, presGPa, color=colour)
        ax.set_ylim((0, 6))
        ax.set_xlim((0, 0.26))
        bx.set_xlim((1180, 1580))
        ax.invert_yaxis()
        for axis in [ax, bx]:
            axis.grid()
            axis.xaxis.tick_top()
        plt.show()
