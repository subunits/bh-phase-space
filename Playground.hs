-- Higher-Dimensional Black Holes -- single-file playground version
-- Paste into https://play.haskell.org (GHC 9.x, base only)
-- Covers: Myers-Perry, Black Rings, Gregory-Laflamme, Phase Space, ZAMOs

module Main where

import Data.List (intercalate)
import Text.Printf (printf)

-- ============================================================
-- CONSTANTS & SPHERE AREAS
-- ============================================================

-- | Omega_n = area of unit n-sphere = 2pi^{(n+1)/2} / Gamma((n+1)/2)
sphereArea :: Int -> Double
sphereArea n = 2 * pi ** (fromIntegral (n+1) / 2) / gammaHalf (n+1)
  where
    -- gammaHalf k = Gamma(k/2)
    -- Base cases: Gamma(1/2) = sqrt(pi), Gamma(1) = 1
    -- Recursion:  Gamma(k/2) = ((k-2)/2) * Gamma((k-2)/2)
    gammaHalf 1 = sqrt pi        -- Gamma(1/2)
    gammaHalf 2 = 1.0            -- Gamma(1)
    gammaHalf k = fromIntegral (k-2) / 2.0 * gammaHalf (k-2)

-- Cusp angular momentum for d=5 ring: j = sqrt(27/32)
jCusp :: Double
jCusp = sqrt (27/32)

-- ============================================================
-- MATH UTILS
-- ============================================================

-- | Bisection root finder. Assumes f(lo)*f(hi) < 0.
bisect :: (Double -> Double) -> Double -> Double -> Double -> Maybe Double
bisect f lo hi tol
  | f lo * f hi > 0 = Nothing
  | otherwise        = Just (go lo hi)
  where
    go a b
      | b - a < tol = (a+b)/2
      | f mid <= 0  = if f a <= 0 then go mid b else go a mid
      | otherwise   = if f a >= 0 then go mid b else go a mid
      where mid = (a+b)/2

linspace :: Double -> Double -> Int -> [Double]
linspace lo hi n = [lo + fromIntegral i * (hi-lo) / fromIntegral (n-1)
                   | i <- [0..n-1]]

-- ============================================================
-- MYERS-PERRY
-- ============================================================

-- | Mass parameter mu = 16piGM / [(d-2) Omega_{d-2}] (G=1)
massParam :: Int -> Double -> Double
massParam d m = 16 * pi * m / (fromIntegral (d-2) * sphereArea (d-2))

-- | Delta(r) = r^2 + a^2 - mu/r^{d-5} (roots = horizons)
deltaMP :: Int -> Double -> Double -> Double -> Double
deltaMP d mu a r = r^2 + a^2 - mu / r ^ (d-5)

-- | d=5 exact horizon: r_H = sqrt(mu - a^2). Nothing if Kerr bound violated.
horizonD5 :: Double -> Double -> Maybe Double
horizonD5 mu a
  | a^2 > mu  = Nothing
  | otherwise = Just (sqrt (mu - a^2))

-- | d>=6 horizon: numerical root of Delta.
horizonDge6 :: Int -> Double -> Double -> Maybe Double
horizonDge6 d mu a =
  bisect (deltaMP d mu a) 1e-9 rMax 1e-10
  where rMax = 10 * (mu ** (1 / fromIntegral (d-3)) + a + 1)

-- | Dispatch horizon by dimension.
horizonMP :: Int -> Double -> Double -> Maybe Double
horizonMP d mu a
  | d == 5    = horizonD5 mu a
  | otherwise = horizonDge6 d mu a

-- | Horizon area: A_H = r_H^{d-4} (r_H^2 + a^2) Omega_{d-2}
areaMP :: Int -> Double -> Double -> Double
areaMP d rH a = rH^(d-4) * (rH^2 + a^2) * sphereArea (d-2)

-- | Parametric a_H(j) curve point at nu = r_H/a
-- Returns (dimensionless area, dimensionless spin)
mpCurvePoint :: Int -> Double -> (Double, Double)
mpCurvePoint d nu =
  ( aHpow ** (1/dm3)
  , jPow  ** (1/dm3) )
  where
    dm3   = fromIntegral (d-3)
    dm4   = fromIntegral (d-4)
    ratio = sphereArea (d-3) / sphereArea (d-2)
    aHpow = 8*pi * (dm4/dm3)**(dm3/2) * ratio * nu^2 / (1+nu^2)
    jPow  = pi / dm3**(dm3/2) * ratio * nu**(5 - fromIntegral d) / (1+nu^2)

-- ============================================================
-- BLACK RING (d=5 only)
-- ============================================================

-- | Equilibrium lambda = 2nu/(1+nu^2) -- removes conical singularity.
eqLambda :: Double -> Double
eqLambda nu = 2*nu / (1 + nu^2)

-- | F(xi) = 1 + lambda*xi
bigF :: Double -> Double -> Double
bigF lam xi = 1 + lam * xi

-- | G(xi) = (1-xi^2)(1+nu*xi)
bigG :: Double -> Double -> Double
bigG nu xi = (1 - xi^2) * (1 + nu*xi)

-- | Ring mass: M = 3pi R^2 lambda / [4G (1-nu)]
ringMass :: Double -> Double -> Double -> Double
ringMass nu lam r = 3*pi*r^2*lam / (4*(1-nu))

-- | Ring angular momentum
ringJ :: Double -> Double -> Double -> Double
ringJ nu lam r =
  pi*r^3 * sqrt (lam*(lam-nu)*(1+lam)) / (2*(1-nu)^2)

-- | Horizon area
ringArea :: Double -> Double -> Double -> Double
ringArea nu lam r =
  8*pi^2*r^3 * nu**(3/2) * sqrt (lam*(1-lam^2))
  / ((1-nu)^2 * (1+nu))

-- | Angular velocity: Omega_H = (1/R) sqrt[(lambda-nu)/(lambda(1+lambda))]
ringOmega :: Double -> Double -> Double -> Double
ringOmega nu lam r = (1/r) * sqrt ((lam-nu) / (lam*(1+lam)))

-- | Parametric a_H(j) curve at nu (equilibrium):
-- a_H = 2sqrt[nu(1-nu)], j = sqrt[(1+nu)^3/8nu]
ringCurvePoint :: Double -> (Double, Double)
ringCurvePoint nu =
  ( 2 * sqrt (nu*(1-nu))
  , sqrt ((1+nu)^3 / (8*nu)) )

-- ============================================================
-- GREGORY-LAFLAMME
-- ============================================================

-- | k_GL = 2pi/r_H -- modes with k < k_GL are unstable.
kGL :: Double -> Double
kGL rH = 2*pi / rH

-- | Is a circular black string of ring-radius R unstable?
glUnstable :: Double -> Double -> Bool
glUnstable bigR rH = bigR > rH

-- | GL growth rate Im(Omega) ~ sqrt(k_GL^2 - k^2) / r_H
glGrowthRate :: Double -> Double -> Maybe Double
glGrowthRate k rH
  | k >= kGL rH = Nothing
  | otherwise   = Just $ sqrt (kGL rH ^ 2 - k^2) / rH

-- ============================================================
-- PHASE SPACE (d=5)
-- ============================================================

-- | Count of distinct d=5 solutions at spin j.
countSolutions :: Double -> Int
countSolutions j
  | j < jCusp = 1   -- MP only
  | j < 1.0   = 3   -- MP + thin + plump ring
  | j == 1.0  = 2
  | otherwise = 1   -- thin ring only

-- ============================================================
-- ZAMO FRAME DRAGGING
-- ============================================================

-- | Omega_ZAMO = -g_{tphi}/g_{phiphi} in equatorial plane, single-spin MP.
zamoOmega :: Int -> Double -> Double -> Double -> Double
zamoOmega d mu a r =
  let rPow   = r ^ (d-5)
      sig    = r^2 + a^2
      gtphi  = -(mu*a) / (rPow * sig)
      gphiphi = (r^2 + a^2) + mu*a^2 / (rPow * sig)
  in -gtphi / gphiphi

-- | 5D equal-spin (a=b) ZAMO from symmetric MP metric:
-- Omega = mu*a / (r^4 + a^2*r^2 + mu*a^2)
-- Valid only when both spin parameters are equal; differs from
-- zamoOmega (single-spin) because the metric components differ.
zamoOmega5DEqual :: Double -> Double -> Double -> Double
zamoOmega5DEqual mu a r = mu*a / (r^4 + a^2*r^2 + mu*a^2)

-- ============================================================
-- MAIN -- demo output
-- ============================================================

main :: IO ()
main = do
  section "n-sphere areas Omega_n"
  mapM_ (\n -> printf "  Omega_%d = %.6f\n" n (sphereArea n)) [0..6::Int]

  section "Myers-Perry d=5 horizons (M=1, G=1)"
  let mu5 = massParam 5 1.0
  printf "  mu = %.6f\n" mu5
  mapM_ (\a -> case horizonD5 mu5 a of
    Nothing -> printf "  a=%.2f  KERR BOUND VIOLATED\n" a
    Just rH -> printf "  a=%.2f  r_H=%.6f  A_H=%.6f\n"
                 a rH (areaMP 5 rH a)
    ) [0.0, 0.3, 0.6, 0.9, 1.0, 1.1]

  section "Myers-Perry d=6 horizons (M=1, G=1) -- no spin bound"
  let mu6 = massParam 6 1.0
  mapM_ (\a -> case horizonDge6 6 mu6 a of
    Nothing -> printf "  a=%-6.0f  no horizon\n" a
    Just rH -> printf "  a=%-6.0f  r_H=%.6f\n" a rH
    ) [0, 1, 10, 100, 1000]

  section "MP a_H(j) curve -- d=5 (j increases toward Kerr bound at j=1)"
  printf "  %10s %10s %s\n" ("j"::String) ("a_H"::String) ("regime"::String)
  let logNus   = map (\t -> 10**(t)) $ reverse $ linspace (log10 1.001) (log10 20.0) 13
      log10 x  = log x / log 10
      labelNu nu = if nu < 1.05 then "near-extremal" else "" :: String
  mapM_ (\nu -> let (aH,j) = mpCurvePoint 5 nu
                in printf "  %10.6f %10.6f %s\n" j aH (labelNu nu)
        ) logNus

  section "MP a_H(j) curve -- d=6 (j unbounded, a_H->0 as spin grows)"
  printf "  %10s %10s\n" ("j"::String) ("a_H"::String)
  mapM_ (\nu -> let (aH,j) = mpCurvePoint 6 nu
                in printf "  %10.6f %10.6f\n" j aH
        ) (reverse $ linspace 0.1 5.0 8)

  section "Black ring -- physical quantities at equilibrium (R=1, G=1)"
  printf "  %-6s %-8s %-10s %-10s %-10s %-8s\n"
    ("nu"::String) ("branch"::String) ("M"::String)
    ("J"::String) ("A_H"::String) ("Omega_H"::String)
  mapM_ (\nu ->
    let lam = eqLambda nu
        m   = ringMass nu lam 1.0
        j   = ringJ    nu lam 1.0
        aH  = ringArea nu lam 1.0
        om  = ringOmega nu lam 1.0
        br  = if nu < 0.5 then "thin" else "plump" :: String
    in printf "  %.2f  %-8s %-10.4f %-10.4f %-10.4f %.4f\n"
         nu br m j aH om
    ) [0.1, 0.3, 0.5, 0.7, 0.9]

  section "Black ring a_H(j) curve -- cusp at nu=0.5, j~0.919"
  printf "  %-6s %-10s %-10s %s\n"
    ("nu"::String) ("j"::String) ("a_H"::String) ("branch"::String)
  mapM_ (\nu ->
    let (aH,j) = ringCurvePoint nu
        br = if nu < 0.5 then "thin" else "plump" :: String
    in printf "  %.3f  %-10.6f %-10.6f %s\n" nu j aH br
    ) (linspace 0.05 0.95 10)

  section "Gregory-Laflamme instability (r_H = 1.0)"
  let rH = 1.0
  printf "  k_GL = %.4f\n" (kGL rH)
  mapM_ (\bigR ->
    printf "  ring R=%.1f  %s\n" bigR
      (if glUnstable bigR rH then "UNSTABLE" else "stable" :: String)
    ) [0.5, 1.0, 1.5, 2.0, 5.0]
  mapM_ (\k -> case glGrowthRate k rH of
    Nothing -> printf "  k=%.2f  stable\n" k
    Just gr -> printf "  k=%.2f  Im(Omega)=%.4f\n" k gr
    ) [1.0, 3.0, 5.0, 6.28, 7.0]

  section "Phase space -- solution count at given j (d=5)"
  mapM_ (\j -> printf "  j=%.3f  %d solution(s)\n" j (countSolutions j))
    [0.5, 0.85, 0.92, 0.95, 0.99, 1.0, 1.1]

  section "ZAMO frame dragging -- d=5 MP, M=1, a=0.5 (single spin)"
  let a = 0.5
  case horizonD5 mu5 a of
    Nothing -> putStrLn "  Kerr bound violated."
    Just rH -> do
      printf "  r_H = %.4f\n" rH
      printf "  %-8s %-20s %-20s\n"
        ("r"::String)
        ("Omega (single-spin)"::String)
        ("Omega (equal-spin a=b)"::String)
      mapM_ (\r ->
        printf "  %-8.3f %-20.6f %-20.6f\n"
          r (zamoOmega 5 mu5 a r) (zamoOmega5DEqual mu5 a r)
        ) [rH*1.01, rH*1.5, rH*2, rH*5, rH*10, 50.0]
      putStrLn "  Note: equal-spin formula assumes a=b; differs from single-spin near horizon."

-- ---------------------------------------------------------------------------
section :: String -> IO ()
section title = do
  putStrLn ""
  putStrLn $ "--- " <> title
  putStrLn $ replicate (length title + 5) '-'
