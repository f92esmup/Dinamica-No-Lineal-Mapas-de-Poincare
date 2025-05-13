import sympy as sp

print("Calculando derivadas simbólicas del Hamiltoniano del péndulo doble...")

# 1. Definir Símbolos (reales para simplificaciones correctas)
theta1, theta2, p1, p2, m, L, g = sp.symbols('theta1 theta2 p1 p2 m L g', real=True, positive=True)
# Asumimos m, L, g positivos

# 2. Definir Expresiones Auxiliares y el Hamiltoniano (H = T + U)
Delta = theta1 - theta2
D_denom = 1 + sp.sin(Delta)**2  # Denominador D = 1 + sin^2(Δ)

# Energía Cinética T
T_num = p1**2 + 2*p2**2 - 2*p1*p2*sp.cos(Delta)
T_den = 2*m*L**2 * D_denom
T = T_num / T_den

# Energía Potencial U
U = m*g*L*(3 - 2*sp.cos(theta1) - sp.cos(theta2))

# Hamiltoniano H
H = T + U
print("\nHamiltoniano H definido simbólicamente.")

# 3. Calcular TODAS las Derivadas Parciales necesarias
print("\nCalculando derivadas parciales...")
dH_dp1 = sp.diff(H, p1)
dH_dp2 = sp.diff(H, p2)
dH_dtheta1 = sp.diff(H, theta1)
dH_dtheta2 = sp.diff(H, theta2)

# 4. Calcular las Derivadas Temporales (Ecuaciones de Hamilton)
theta1_dot_expr = dH_dp1
theta2_dot_expr = dH_dp2
p1_dot_expr = -dH_dtheta1
p2_dot_expr = -dH_dtheta2

# 5. Imprimir resultados
print("-" * 40)
print("Expresión simbólica para theta1_dot = dH/dp1:")
sp.pprint(theta1_dot_expr, wrap_line=False)
print("-" * 40)
print("\nExpresión simbólica para theta2_dot = dH/dp2:")
sp.pprint(theta2_dot_expr, wrap_line=False)
print("-" * 40)

print("\nExpresión simbólica para p1_dot = -dH/dtheta1:")
sp.pprint(p1_dot_expr, wrap_line=False)
print("-" * 40)

print("\nExpresión simbólica para p2_dot = -dH/dtheta2:")
sp.pprint(p2_dot_expr, wrap_line=False)
print("-" * 40)
print("\nDerivadas parciales calculadas.")