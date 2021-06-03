#!/usr/bin/env python3

def func(X):
    from numpy import array, clip, exp, sqrt
    from numpy.linalg import norm
    from scipy.optimize import root
    from scipy.special import erf
    from Constants import (conDepth, d, domainDim, k, kappa, mantTemp, oceAge,
                           pacman, rhoH0, stepLoc, stepWidth, surfTemp)

    def solvStep(var, dist2shall, Z, width, surfTemp, mantTemp, kappa, oceAge,
                 k1, k2, rhoH0, d):
        T, zCON, zOCE = var
        return ((erf(dist2shall / width * 4. - 2.) + 1.) / 2.
                * (zCON - zOCE) + zOCE - Z,
                T - surfTemp - (mantTemp - surfTemp)
                * erf(zOCE / 2. / sqrt(kappa * oceAge)),
                T - (k1 * zCON + k2 - rhoH0 * d ** 2. * exp(-zCON / d)) / 3.)

    depth = clip(domainDim[2] - X[2], 0., domainDim[2])
    k2 = k * surfTemp + rhoH0 * d ** 2
    k1 = (k * mantTemp + rhoH0 * d ** 2 * exp(-conDepth / d) - k2) / conDepth
    conTemp = (k1 * depth + k2 - rhoH0 * d ** 2. * exp(-depth / d)) / k
    oceTemp = surfTemp + (mantTemp - surfTemp) * erf(depth / 2.
                                                     / sqrt(kappa * oceAge))
    # Depth
    if depth > conDepth:
        return mantTemp
    # Ocean
    elif (X[0] <= stepLoc[0][0] - stepWidth / 2
          or X[0] >= stepLoc[0][1] + stepWidth / 2):
        return oceTemp
    elif (X[1] <= stepLoc[1][0] - stepWidth / 2
          or X[1] >= stepLoc[1][1] + stepWidth / 2):
        return oceTemp
    # Outside Steps
    elif (X[1] < stepLoc[1][0] + stepWidth / 2
          and X[0] >= stepLoc[0][0] + stepWidth / 2
          and X[0] <= stepLoc[0][1] - stepWidth / 2):
        res = root(solvStep, x0=(900., 8e4, 4e4),
                   args=(X[1] - stepLoc[1][0] + stepWidth / 2, depth,
                         stepWidth, surfTemp, mantTemp, kappa, oceAge, k1, k2,
                         rhoH0, d))
        return res.x[0]
    elif (X[1] > stepLoc[1][1] - stepWidth / 2
          and X[0] >= stepLoc[0][0] + stepWidth / 2
          and X[0] <= stepLoc[0][1] - stepWidth / 2):
        res = root(solvStep, x0=(900., 8e4, 4e4),
                   args=(stepLoc[1][1] + stepWidth / 2 - X[1], depth,
                         stepWidth, surfTemp, mantTemp, kappa, oceAge, k1, k2,
                         rhoH0, d))
        return res.x[0]
    elif (X[0] < stepLoc[0][0] + stepWidth / 2
          and X[1] >= stepLoc[1][0] + stepWidth / 2
          and X[1] <= stepLoc[1][1] - stepWidth / 2):
        res = root(solvStep, x0=(900., 8e4, 4e4),
                   args=(X[0] - stepLoc[0][0] + stepWidth / 2, depth,
                         stepWidth, surfTemp, mantTemp, kappa, oceAge, k1, k2,
                         rhoH0, d))
        return res.x[0]
    elif (X[0] > stepLoc[0][1] - stepWidth / 2
          and ((X[1] >= stepLoc[1][0] + stepWidth / 2
                and X[1] <= stepLoc[1][0] + (stepLoc[1][1] - stepLoc[1][0]
                                             - pacman - stepWidth) / 2)
               or (X[1] <= stepLoc[1][1] - stepWidth / 2
                   and X[1] >= stepLoc[1][1] - (stepLoc[1][1] - stepLoc[1][0]
                                                - pacman - stepWidth) / 2))):
        res = root(solvStep, x0=(900., 8e4, 4e4),
                   args=(stepLoc[0][1] + stepWidth / 2 - X[0], depth,
                         stepWidth, surfTemp, mantTemp, kappa, oceAge, k1, k2,
                         rhoH0, d))
        return res.x[0]
    # Outside Corners
    elif (X[0] < stepLoc[0][0] + stepWidth / 2
          and X[1] < stepLoc[1][0] + stepWidth / 2):
        dist = norm(array(X[:-1]) - array([stepLoc[0][0] + stepWidth / 2,
                                           stepLoc[1][0] + stepWidth / 2]))
        if dist > stepWidth:
            return oceTemp
        else:
            res = root(solvStep, x0=(900., 8e4, 4e4),
                       args=(stepWidth - dist, depth,
                             stepWidth, surfTemp, mantTemp, kappa, oceAge,
                             k1, k2, rhoH0, d))
            return res.x[0]
    elif (X[0] < stepLoc[0][0] + stepWidth / 2
          and X[1] > stepLoc[1][1] - stepWidth / 2):
        dist = norm(array(X[:-1]) - array([stepLoc[0][0] + stepWidth / 2,
                                           stepLoc[1][1] - stepWidth / 2]))
        if dist > stepWidth:
            return oceTemp
        else:
            res = root(solvStep, x0=(900., 8e4, 4e4),
                       args=(stepWidth - dist, depth,
                             stepWidth, surfTemp, mantTemp, kappa, oceAge,
                             k1, k2, rhoH0, d))
            return res.x[0]
    elif (X[0] > stepLoc[0][1] - stepWidth / 2
          and X[1] < stepLoc[1][0] + stepWidth / 2):
        dist = norm(array(X[:-1]) - array([stepLoc[0][1] - stepWidth / 2,
                                           stepLoc[1][0] + stepWidth / 2]))
        if dist > stepWidth:
            return oceTemp
        else:
            res = root(solvStep, x0=(900., 8e4, 4e4),
                       args=(stepWidth - dist, depth,
                             stepWidth, surfTemp, mantTemp, kappa, oceAge,
                             k1, k2, rhoH0, d))
            return res.x[0]
    elif (X[0] > stepLoc[0][1] - stepWidth / 2
          and X[1] > stepLoc[1][1] - stepWidth / 2):
        dist = norm(array(X[:-1]) - array([stepLoc[0][1] - stepWidth / 2,
                                           stepLoc[1][1] - stepWidth / 2]))
        if dist > stepWidth:
            return oceTemp
        else:
            res = root(solvStep, x0=(900., 8e4, 4e4),
                       args=(stepWidth - dist, depth,
                             stepWidth, surfTemp, mantTemp, kappa, oceAge,
                             k1, k2, rhoH0, d))
            return res.x[0]
    # Inside Steps
    elif (X[0] > stepLoc[0][1] - pacman - stepWidth / 2
          and X[0] < stepLoc[0][1] - pacman + stepWidth / 2
          and X[1] <= stepLoc[1][1] - (stepLoc[1][1] - stepLoc[1][0]
                                       - pacman + stepWidth) / 2
          and X[1] >= stepLoc[1][0] + (stepLoc[1][1] - stepLoc[1][0]
                                       - pacman + stepWidth) / 2):
        res = root(solvStep, x0=(900., 8e4, 4e4),
                   args=(stepWidth - X[0] + stepLoc[0][1] - pacman
                         - stepWidth / 2, depth, stepWidth, surfTemp, mantTemp,
                         kappa, oceAge, k1, k2, rhoH0, d))
        return res.x[0]
    elif (X[0] >= stepLoc[0][1] - pacman + stepWidth / 2
          and X[0] <= stepLoc[0][1] - stepWidth / 2
          and X[1] < stepLoc[1][0] + (stepLoc[1][1] - stepLoc[1][0]
                                      - pacman + stepWidth) / 2
          and X[1] > stepLoc[1][0] + (stepLoc[1][1] - stepLoc[1][0]
                                      - pacman - stepWidth) / 2):
        res = root(solvStep, x0=(900., 8e4, 4e4),
                   args=(stepWidth - X[1] + stepLoc[1][0]
                         + (stepLoc[1][1] - stepLoc[1][0] - pacman
                            - stepWidth) / 2, depth, stepWidth, surfTemp,
                         mantTemp, kappa, oceAge, k1, k2, rhoH0, d))
        return res.x[0]
    elif (X[0] >= stepLoc[0][1] - pacman + stepWidth / 2
          and X[0] <= stepLoc[0][1] - stepWidth / 2
          and X[1] < stepLoc[1][1] - (stepLoc[1][1] - stepLoc[1][0]
                                      - pacman - stepWidth) / 2
          and X[1] > stepLoc[1][1] - (stepLoc[1][1] - stepLoc[1][0]
                                      - pacman + stepWidth) / 2):
        res = root(solvStep, x0=(900., 8e4, 4e4),
                   args=(X[1] - stepLoc[1][1] + (stepLoc[1][1] - stepLoc[1][0]
                                                 - pacman + stepWidth) / 2,
                         depth, stepWidth, surfTemp, mantTemp, kappa, oceAge,
                         k1, k2, rhoH0, d))
        return res.x[0]
    # Inner Corners
    elif (X[0] > stepLoc[0][1] - stepWidth / 2
          and X[1] > stepLoc[1][1] - (stepLoc[1][1] - stepLoc[1][0]
                                      - pacman + stepWidth) / 2):
        dist = norm(array(X[:-1]) - array([stepLoc[0][1] - stepWidth / 2,
                                           stepLoc[1][1]
                                           - (stepLoc[1][1] - stepLoc[1][0]
                                              - pacman - stepWidth) / 2]))
        if dist > stepWidth:
            return oceTemp
        else:
            res = root(solvStep, x0=(900., 8e4, 4e4),
                       args=(stepWidth - dist, depth,
                             stepWidth, surfTemp, mantTemp, kappa, oceAge,
                             k1, k2, rhoH0, d))
            return res.x[0]
    elif (X[0] > stepLoc[0][1] - stepWidth / 2
          and X[1] < stepLoc[1][0] + (stepLoc[1][1] - stepLoc[1][0]
                                      - pacman + stepWidth) / 2):
        dist = norm(array(X[:-1]) - array([stepLoc[0][1] - stepWidth / 2,
                                           stepLoc[1][0]
                                           + (stepLoc[1][1] - stepLoc[1][0]
                                              - pacman - stepWidth) / 2]))
        if dist > stepWidth:
            return oceTemp
        else:
            res = root(solvStep, x0=(900., 8e4, 4e4),
                       args=(stepWidth - dist, depth,
                             stepWidth, surfTemp, mantTemp, kappa, oceAge,
                             k1, k2, rhoH0, d))
            return res.x[0]
    elif (X[0] > stepLoc[0][1] - pacman - stepWidth / 2
          and X[1] < stepLoc[1][1] - (stepLoc[1][1] - stepLoc[1][0]
                                      - pacman - stepWidth) / 2
          and X[1] > stepLoc[1][1] - (stepLoc[1][1] - stepLoc[1][0]
                                      - pacman + stepWidth) / 2):
        dist = norm(array(X[:-1]) - array([stepLoc[0][1] - pacman
                                           + stepWidth / 2,
                                           stepLoc[1][1]
                                           - (stepLoc[1][1] - stepLoc[1][0]
                                              - pacman + stepWidth) / 2]))
        if dist > stepWidth:
            return conTemp
        else:
            res = root(solvStep, x0=(900., 8e4, 4e4),
                       args=(dist, depth,
                             stepWidth, surfTemp, mantTemp, kappa, oceAge,
                             k1, k2, rhoH0, d))
            return res.x[0]
    elif (X[0] > stepLoc[0][1] - pacman - stepWidth / 2
          and X[1] < stepLoc[1][0] + (stepLoc[1][1] - stepLoc[1][0]
                                      - pacman + stepWidth) / 2
          and X[1] > stepLoc[1][0] + (stepLoc[1][1] - stepLoc[1][0]
                                      - pacman - stepWidth) / 2):
        dist = norm(array(X[:-1]) - array([stepLoc[0][1] - pacman
                                           + stepWidth / 2,
                                           stepLoc[1][0]
                                           + (stepLoc[1][1] - stepLoc[1][0]
                                              - pacman + stepWidth) / 2]))
        if dist > stepWidth:
            return conTemp
        else:
            res = root(solvStep, x0=(900., 8e4, 4e4),
                       args=(dist, depth,
                             stepWidth, surfTemp, mantTemp, kappa, oceAge,
                             k1, k2, rhoH0, d))
            return res.x[0]
    # Ocean Indent
    elif (X[0] >= stepLoc[0][1] - pacman + stepWidth / 2
          and X[1] <= stepLoc[1][1] - (stepLoc[1][1] - stepLoc[1][0]
                                       - pacman + stepWidth) / 2
          and X[1] >= stepLoc[1][0] + (stepLoc[1][1] - stepLoc[1][0]
                                       - pacman + stepWidth) / 2):
        return oceTemp
    # Continent
    else:
        return conTemp
