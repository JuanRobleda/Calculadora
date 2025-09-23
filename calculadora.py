# app.py
import streamlit as st
from fractions import Fraction
import base64
import os
import matplotlib.pyplot as plt

# -------------------------
# Helpers: Fracciones y LaTeX
# -------------------------
def frac_str(q: Fraction) -> str:
    if isinstance(q, Fraction):
        if q.denominator == 1:
            return str(q.numerator)
        return f"{q.numerator}/{q.denominator}"
    return str(q)

def frac_to_latex(f):
    if isinstance(f, Fraction):
        if f.denominator == 1:
            return str(f.numerator)
        else:
            return r"\tfrac{" + str(f.numerator) + "}{" + str(f.denominator) + "}"
    return str(f)

def latex_matrix(M):
    rows = []
    for row in M:
        elems = [frac_to_latex(x) for x in row]
        rows.append(" & ".join(elems))
    body = r" \\\\ ".join(rows)
    return r"\begin{bmatrix}" + body + r"\end{bmatrix}"

def latex_column_vector(v):
    rows = [frac_to_latex(x) for x in v]
    body = r" \\\\ ".join(rows)
    return r"\begin{bmatrix}" + body + r"\end{bmatrix}"

def latex_matrix_highlight(M, pivotes=[]):
    rows = []
    for i, row in enumerate(M):
        elems = []
        for j, val in enumerate(row):
            latex_val = frac_to_latex(val)
            if (i, j) in pivotes:
                elems.append(r"\color{red}{" + latex_val + "}")
            else:
                elems.append(latex_val)
        rows.append(" & ".join(elems))
    body = r" \\\\ ".join(rows)
    return r"\begin{bmatrix}" + body + r"\end{bmatrix}"

def pretty_equations(A_rows, B_parts, var_symbol="x"):
    """
    Construye lista de strings LaTeX que representan el sistema de ecuaciones
    derivado de A_rows * x = B_parts.
    A_rows: lista de filas (cada fila lista de Fraction)
    B_parts: lista de Fraction (longitud igual a filas)
    Devuelve lista de cadenas LaTeX.
    """
    eqs = []
    n = len(A_rows)
    if n == 0:
        return eqs
    m = len(A_rows[0])
    for i in range(n):
        terms = []
        for j in range(m):
            a = A_rows[i][j]
            if a == 0:
                continue
            coeff = frac_to_latex(a)
            # handle sign and coefficient formatting
            # if coefficient is 1 or -1, show sign and variable accordingly
            if isinstance(a, Fraction) and a == 1:
                term = f"{var_symbol}_{{{j+1}}}"
            elif isinstance(a, Fraction) and a == -1:
                term = f"-{var_symbol}_{{{j+1}}}"
            else:
                term = f"{coeff}{var_symbol}_{{{j+1}}}"
            terms.append(term)
        if not terms:
            left = "0"
        else:
            # join terms with proper +/-
            # terms already include minus if coefficient negative
            left = " + ".join(terms)
            # fix patterns like "+ -" -> "-"
            left = left.replace("+ -", "- ")
        right = frac_to_latex(B_parts[i])
        eqs.append(left + " = " + right)
    return eqs

# -------------------------
# GAUSS y GAUSS-JORDAN
# -------------------------
def gauss_jordan_steps(M):
    """
    Espera M: lista de filas (cada fila: lista de Fraction o numéricos).
    Devuelve lista de pasos: (descripcion, matriz_actual (deep copy de Fractions)).
    Incluye la matriz inicial.
    """
    n = len(M)
    m = len(M[0]) if n > 0 else 0
    A = [[x if isinstance(x, Fraction) else Fraction(x) for x in fila] for fila in M]
    pasos = []
    pasos.append(("Matriz inicial", [row.copy() for row in A]))
    r = 0
    for c in range(m - 1):
        if r >= n:
            break
        piv = next((i for i in range(r, n) if A[i][c] != 0), None)
        if piv is None:
            continue
        if piv != r:
            A[r], A[piv] = A[piv], A[r]
            pasos.append((f"Intercambio F{r+1} ↔ F{piv+1}", [row.copy() for row in A]))
        k = A[r][c]
        if k != 1 and k != 0:
            A[r] = [x / k for x in A[r]]
            pasos.append((f"F{r+1} → F{r+1}/{frac_str(k)}", [row.copy() for row in A]))
        for i in range(n):
            if i == r:
                continue
            factor = A[i][c]
            if factor != 0:
                A[i] = [A[i][j] - factor * A[r][j] for j in range(m)]
                pasos.append((f"F{i+1} → F{i+1} - ({frac_str(factor)})·F{r+1}", [row.copy() for row in A]))
        r += 1
    return pasos

def resolver_sistema(M, return_steps=False):
    pasos = gauss_jordan_steps(M)
    if not pasos:
        if return_steps:
            return pasos, "⚠️ No se pudo reducir la matriz (o ya está en RREF).", [], []
        return "⚠️ No se pudo reducir la matriz (o ya está en RREF).", [], []
    R = pasos[-1][1]
    texto = "✅ Matriz final (RREF):\n" + str(R)
    n, m = len(R), len(R[0])
    pivot_cols = []
    for i in range(n):
        for j in range(m - 1):
            if R[i][j] != 0:
                pivot_cols.append(j)
                break
    pivot_cols = sorted(set(pivot_cols))
    free_vars = [j for j in range(m - 1) if j not in pivot_cols]
    texto += "\n\n📌 Variables básicas: " + (", ".join([f"x{j+1}" for j in pivot_cols]) if pivot_cols else "No hay variables básicas.")
    texto += "\n📌 Variables libres: " + (", ".join([f"x{j+1}" for j in free_vars]) if free_vars else "No hay variables libres.")
    if free_vars:
        texto += "\n\n⚡ Sistema con solución general parametrizada:\n"
        soluciones = ["0"] * (m - 1)
        for i in range(n):
            piv_col = next((j for j in range(m - 1) if R[i][j] != 0), None)
            if piv_col is not None:
                expr = []
                for j in free_vars:
                    if R[i][j] != 0:
                        expr.append(f"-({frac_str(R[i][j])})*x{j+1}")
                rhs = frac_str(R[i][-1])
                soluciones[piv_col] = rhs + (" + " + " + ".join(expr) if expr else "")
        for j in free_vars:
            soluciones[j] = f"x{j+1}"
        for idx, val in enumerate(soluciones):
            texto += f"x{idx+1} = {val}\n"
    else:
        texto += "\n\n✅ Sistema consistente con única solución:\n"
        sol = []
        for i in range(n):
            piv_col = next((j for j in range(m - 1) if R[i][j] != 0), None)
            if piv_col is not None:
                sol.append((piv_col, R[i][-1]))
        for idx, val in sorted(sol):
            texto += f"x{idx+1} = {frac_str(val)}\n"
    if return_steps:
        return pasos, texto, pivot_cols, free_vars
    return texto, pivot_cols, free_vars

# -------------------------
# Fondo en video
# -------------------------
def set_video_background(video_file):
    if not os.path.exists(video_file):
        return
    with open(video_file, "rb") as f:
        data = f.read()
    b64 = base64.b64encode(data).decode()
    video_html = f"""
    <style>
    .stApp {{
        background: transparent;
        font-family: 'Arial', sans-serif;
    }}
    .resultado {{
        background-color: white;
        color: black;
        padding: 15px;
        border-radius: 10px;
        margin-top: 15px;
    }}
    </style>
    <video autoplay muted loop id="bgvideo"
        style="position: fixed; right: 0; bottom: 0; min-width: 100%; min-height: 100%;
        z-index: -1; opacity: 0.2;">
        <source src="data:video/mp4;base64,{b64}" type="video/mp4">
    </video>
    """
    st.markdown(video_html, unsafe_allow_html=True)

# -------------------------
# Operaciones con vectores
# -------------------------
def parse_vector_text(s):
    parts = [p.strip() for p in s.replace(",", " ").split() if p.strip() != ""]
    vec = []
    for p in parts:
        try:
            vec.append(Fraction(p))
        except Exception:
            vec.append(p)
    return vec

def suma_vectores_varios(vectors):
    n = len(vectors[0])
    resultado = []
    for i in range(n):
        comp_vals = [v[i] for v in vectors]
        if all(isinstance(x, Fraction) for x in comp_vals):
            s = sum(comp_vals, Fraction(0))
            resultado.append(s)
        else:
            resultado.append(" + ".join(frac_to_latex(x) for x in comp_vals))
    return resultado

def resta_vectores(v1, v2):
    return [a - b for a, b in zip(v1, v2)]

def producto_punto(v1, v2):
    if all(isinstance(a, Fraction) and isinstance(b, Fraction) for a, b in zip(v1, v2)):
        return sum(a * b for a, b in zip(v1, v2))
    productos = [f"{frac_to_latex(a)}\\cdot{frac_to_latex(b)}" for a, b in zip(v1, v2)]
    return " + ".join(productos)

# -------------------------
# APP
# -------------------------
st.set_page_config(page_title="Calculadora Álgebra Lineal", page_icon="🧮", layout="wide")
set_video_background("fondo.mp4")

st.title(" Calculadora de Álgebra Lineal ")
st.write("✨ Usa la barra lateral para elegir qué operación deseas realizar.")

seccion = st.sidebar.selectbox("Selecciona sección", [
    "Vectores (varios)",
    "Combinación lineal (AX)",
    "Ecuación matricial AX = B",
    "Gauss / Gauss-Jordan (RREF)"
])

# -------------------------
# SECCIÓN: VECTORES (varios)
# -------------------------
if seccion == "Vectores (varios)":
    st.header("📐 Operaciones con Vectores")
    n_vectores = st.number_input("¿Cuántos vectores deseas ingresar?", min_value=2, max_value=10, value=3, step=1)
    vectores = []
    for i in range(n_vectores):
        texto = st.text_input(f"Vector {i+1}", value="1 2", key=f"vec_{i}")
        vectores.append(parse_vector_text(texto))

    if st.button("Ver resultado"):
        dimensiones = [len(v) for v in vectores]
        if len(set(dimensiones)) != 1:
            st.error("⚠️ Todos los vectores deben tener la misma dimensión.")
        else:
            st.markdown('<div class="resultado">', unsafe_allow_html=True)
            for i, v in enumerate(vectores, 1):
                st.latex(r"\vec{v_" + str(i) + "} = " + latex_column_vector(v))

            st.write("➕ **Suma de todos los vectores (detalle por componente):**")
            pasos_componentes = []
            suma_res = []
            for idx in range(len(vectores[0])):
                comp_vals = [v[idx] for v in vectores]
                sumandos = " + ".join(frac_to_latex(x) for x in comp_vals)
                if all(isinstance(x, Fraction) for x in comp_vals):
                    suma_val = sum(comp_vals, Fraction(0))
                    pasos_componentes.append(f"{sumandos} = {frac_to_latex(suma_val)}")
                    suma_res.append(suma_val)
                else:
                    pasos_componentes.append(f"{sumandos}")
                    suma_res.append(sumandos)
            pasos_latex = r"\begin{bmatrix}" + r" \\\\ ".join(pasos_componentes) + r"\end{bmatrix}"
            st.latex(r"\vec{S} = " + pasos_latex)
            st.subheader("Resultado final de la suma")
            st.latex(r"\vec{S} = " + latex_column_vector(suma_res))

            st.write("Detalle: suma componente a componente:")
            for idx, texto_comp in enumerate(pasos_componentes, start=1):
                st.latex(r"x_{" + str(idx) + r"}: \quad " + texto_comp)

            st.write("➖ **Resta por pares (detalle):**")
            for i in range(len(vectores)-1):
                v1, v2 = vectores[i], vectores[i+1]
                pasos_r = []
                for a, b in zip(v1, v2):
                    try:
                        pasos_r.append(frac_to_latex(a) + " - " + frac_to_latex(b) + " = " + frac_to_latex(a-b))
                    except Exception:
                        pasos_r.append(frac_to_latex(a) + " - " + frac_to_latex(b))
                pasos_latex_r = r"\begin{bmatrix}" + r" \\\\ ".join(pasos_r) + r"\end{bmatrix}"
                st.latex(r"\vec{v_" + str(i+1) + "} - \vec{v_" + str(i+2) + "} = " + pasos_latex_r)
                try:
                    resultado_resta = [a - b for a, b in zip(v1, v2)]
                    st.latex(r"= " + latex_column_vector(resultado_resta))
                except Exception:
                    pass

            st.write("✖ **Producto punto (pares) con detalle por componentes:**")
            for i in range(len(vectores)-1):
                v1, v2 = vectores[i], vectores[i+1]
                productos = []
                tiene_todos_numericos = all(isinstance(a, Fraction) and isinstance(b, Fraction) for a, b in zip(v1, v2))
                for a, b in zip(v1, v2):
                    if isinstance(a, Fraction) and isinstance(b, Fraction):
                        productos.append(f"{frac_to_latex(a)} \\cdot {frac_to_latex(b)} = {frac_to_latex(a*b)}")
                    else:
                        productos.append(f"{frac_to_latex(a)} \\cdot {frac_to_latex(b)}")
                productos_latex = r"\begin{bmatrix}" + r" \\\\ ".join(productos) + r"\end{bmatrix}"
                st.latex(r"\text{Componentes: } " + productos_latex)
                if tiene_todos_numericos:
                    total = sum(a*b for a, b in zip(v1, v2))
                    suma_productos = " + ".join(frac_to_latex(a*b) for a, b in zip(v1, v2))
                    st.latex(r"\vec{v_" + str(i+1) + "} \cdot \vec{v_" + str(i+2) + "} = " + suma_productos + " = " + frac_to_latex(total))
                else:
                    prod_expr = " + ".join(frac_to_latex(a) + r"\cdot" + frac_to_latex(b) for a, b in zip(v1, v2))
                    st.latex(r"\vec{v_" + str(i+1) + "} \cdot \vec{v_" + str(i+2) + "} = " + prod_expr)

            st.markdown('</div>', unsafe_allow_html=True)

            if len(vectores[0]) == 2 and all(all(isinstance(x, Fraction) for x in v) for v in vectores):
                colores = ["blue", "red", "green", "orange", "purple", "brown", "pink", "gray", "cyan", "black"]
                fig, ax = plt.subplots()
                ax.axhline(0, color='gray', linewidth=1)
                ax.axvline(0, color='gray', linewidth=1)
                ax.set_aspect('equal', adjustable='box')
                ax.grid(True, linestyle='--', alpha=0.5)
                for i, v in enumerate(vectores):
                    ax.quiver(0, 0, float(v[0]), float(v[1]),
                              angles='xy', scale_units='xy', scale=1,
                              color=colores[i % len(colores)], label=f"v{i+1}")
                max_x = max(abs(float(v[0])) for v in vectores) + 1
                max_y = max(abs(float(v[1])) for v in vectores) + 1
                ax.set_xlim(-max_x, max_x)
                ax.set_ylim(-max_y, max_y)
                ax.legend()
                ax.set_xlabel("$x_1$")
                ax.set_ylabel("$x_2$")
                st.subheader("📊 Gráfico de vectores individuales")
                st.pyplot(fig)

                suma_parcial = [Fraction(0), Fraction(0)]
                fig2, ax2 = plt.subplots()
                ax2.axhline(0, color='gray', linewidth=1)
                ax2.axvline(0, color='gray', linewidth=1)
                ax2.set_aspect('equal', adjustable='box')
                ax2.grid(True, linestyle='--', alpha=0.5)
                for i, v in enumerate(vectores):
                    ax2.quiver(float(suma_parcial[0]), float(suma_parcial[1]),
                               float(v[0]), float(v[1]),
                               angles='xy', scale_units='xy', scale=1,
                               color=colores[i % len(colores)], label=f"Paso {i+1}: +v{i+1}")
                    suma_parcial[0] += v[0]
                    suma_parcial[1] += v[1]
                ax2.set_xlim(-max_x, max_x)
                ax2.set_ylim(-max_y, max_y)
                ax2.legend()
                ax2.set_xlabel("$x_1$")
                ax2.set_ylabel("$x_2$")
                st.subheader("📊 Gráfico de la suma paso a paso")
                st.pyplot(fig2)

# -------------------------
# SECCIÓN: COMBINACIÓN LINEAL (AX)
# -------------------------
elif seccion == "Combinación lineal (AX)":
    st.header("🔄 Combinación lineal AX")
    st.write("Tienes dos modos en esta sección:")
    st.write("- **Calcular A·X** (entrada X simbólica o numérica).")
    st.write("- **Verificar si b es combinación lineal** de las columnas de A (muestra procedimiento Gauss-Jordan).")

    modo = st.radio("Selecciona modo:", ["Calcular A·X", "Verificar combinación lineal"])
    A_text = st.text_area("Matriz A (filas separadas por saltos de línea)", value="1 0\n-2 5\n-5 6", height=120)

    if modo == "Calcular A·X":
        x_text = st.text_input("Vector X (ej: x1 x2)", value="x1 x2")
        if st.button("Ver resultado - A·X"):
            A_rows = []
            for row in A_text.splitlines():
                if row.strip():
                    nums = [p.strip() for p in row.replace(",", " ").split() if p.strip()]
                    fila = []
                    for p in nums:
                        try:
                            fila.append(Fraction(p))
                        except Exception:
                            fila.append(p)
                    A_rows.append(fila)
            x_vec = parse_vector_text(x_text)
            if not A_rows:
                st.error("Introduce una matriz A válida.")
                st.stop()
            if any(len(r) != len(A_rows[0]) for r in A_rows):
                st.error("Todas las filas de A deben tener la misma cantidad de elementos.")
                st.stop()
            if len(x_vec) != len(A_rows[0]):
                st.warning(f"A tiene {len(A_rows[0])} columnas y X tiene {len(x_vec)} elementos. Se hará ajuste automático.")
                if len(x_vec) < len(A_rows[0]):
                    x_vec += [Fraction(0)] * (len(A_rows[0]) - len(x_vec))
                else:
                    x_vec = x_vec[:len(A_rows[0])]
            resultado = []
            pasos_detalle = []
            for i, fila in enumerate(A_rows, start=1):
                expr = []
                suma_val = Fraction(0)
                symbolic = False
                for a, x in zip(fila, x_vec):
                    if isinstance(a, Fraction) and isinstance(x, Fraction):
                        expr.append(f"{frac_to_latex(a)} \\cdot {frac_to_latex(x)}")
                        suma_val += a * x
                    else:
                        expr.append(f"{frac_to_latex(a)} \\cdot {x}")
                        symbolic = True
                if symbolic:
                    pasos_detalle.append(" + ".join(expr))
                    resultado.append(" + ".join(expr))
                else:
                    pasos_detalle.append(" + ".join(expr) + f" = {frac_to_latex(suma_val)}")
                    resultado.append(suma_val)
            st.markdown('<div class="resultado">', unsafe_allow_html=True)
            st.latex(r"A = " + latex_matrix(A_rows))
            st.latex(r"X = " + latex_column_vector(x_vec))
            st.subheader("📝 Pasos del cálculo (fila por fila)")
            pasos_latex = r"\begin{bmatrix}" + r" \\\\ ".join(pasos_detalle) + r"\end{bmatrix}"
            st.latex(r"A \cdot X = " + pasos_latex)
            st.subheader("✅ Resultado final")
            st.latex(r"A \cdot X = " + latex_column_vector(resultado))
            st.markdown('</div>', unsafe_allow_html=True)

    else:  # Verificar combinación lineal
        B_text = st.text_input("Vector B (ej: 7 4 -3) — se verifica si b ∈ span{columnas de A}", value="7 4 -3")
        if st.button("Verificar combinación lineal"):
            A_rows = []
            try:
                for row in A_text.splitlines():
                    if row.strip():
                        parts = [p.strip() for p in row.replace(",", " ").split() if p.strip()]
                        fila = [Fraction(p) for p in parts]
                        A_rows.append(fila)
                B_parts = [Fraction(p) for p in B_text.replace(",", " ").split() if p.strip()]
            except Exception:
                st.error("⚠️ A y B deben contener números enteros o fracciones (ej. 3, -1/2).")
                st.stop()
            if not A_rows:
                st.error("Introduce una matriz A válida.")
                st.stop()
            if any(len(row) != len(A_rows[0]) for row in A_rows):
                st.error("⚠️ Todas las filas de A deben tener la misma longitud.")
                st.stop()
            n_filas_A = len(A_rows)
            n_cols_A = len(A_rows[0])
            if len(B_parts) < n_filas_A:
                B_parts = B_parts + [Fraction(0)] * (n_filas_A - len(B_parts))
            elif len(B_parts) > n_filas_A:
                st.warning("B tiene más entradas que filas de A: se truncará B a las primeras filas.")
                B_parts = B_parts[:n_filas_A]

            M_aug = [row + [B_parts[i]] for i, row in enumerate(A_rows)]

            st.markdown('<div class="resultado">', unsafe_allow_html=True)
            st.latex(r"A = " + latex_matrix(A_rows))
            st.latex(r"b = " + latex_column_vector(B_parts))

            # Mostrar sistema de ecuaciones en notación clásica (como en tu foto)
            st.subheader("📑 Sistema de ecuaciones lineales (forma expandida)")
            eqs = pretty_equations(A_rows, B_parts, var_symbol="x")
            for e in eqs:
                st.latex(e)

            st.subheader("📋 Matriz aumentada [A | b] y procedimiento (Gauss-Jordan)")
            pasos = gauss_jordan_steps(M_aug)
            if pasos:
                for i, (op, mat) in enumerate(pasos, start=1):
                    st.markdown(f"**Paso {i}: {op}**")
                    st.latex(latex_matrix_highlight(mat))
                R = pasos[-1][1]
            else:
                R = [row.copy() for row in M_aug]
                st.write("No se realizaron pasos (la matriz pudo estar ya en RREF).")
                st.latex(latex_matrix(R))

            st.subheader("✅ Matriz final (RREF):")
            st.latex(latex_matrix(R))

            inconsistent = False
            for row in R:
                left = row[:-1]
                right = row[-1]
                if all(x == 0 for x in left) and right != 0:
                    inconsistent = True
                    break

            if inconsistent:
                st.error("❌ El sistema es inconsistente — b NO es combinación lineal de las columnas de A.")
            else:
                n = len(R)
                m = len(R[0])
                pivot_cols = []
                for i in range(n):
                    for j in range(m-1):
                        if R[i][j] != 0:
                            pivot_cols.append(j)
                            break
                pivot_cols = sorted(set(pivot_cols))
                free_vars = [j for j in range(m-1) if j not in pivot_cols]

                if not free_vars:
                    sol = [None] * (m-1)
                    pivot_row_for_col = {}
                    for i in range(n):
                        for j in range(m-1):
                            if R[i][j] != 0:
                                pivot_row_for_col[j] = i
                                break
                    for j in pivot_cols:
                        sol[j] = R[pivot_row_for_col[j]][-1]
                    st.success("✅ b SÍ es combinación lineal de las columnas de A (solución única).")
                    st.subheader("📌 Coeficientes x encontrados:")
                    for idx, val in enumerate(sol):
                        st.write(f"x{idx+1} = {frac_str(val)}")

                    columnas_A = []
                    for col_idx in range(n_cols_A):
                        col = [A_rows[row_idx][col_idx] for row_idx in range(n_filas_A)]
                        columnas_A.append(col)
                    b_reconstruido = [Fraction(0)] * n_filas_A
                    for j, coef in enumerate(sol):
                        if coef is not None:
                            for i in range(n_filas_A):
                                b_reconstruido[i] += coef * columnas_A[j][i]

                    st.markdown("**Combinación lineal encontrada (notación):**")
                    terms = []
                    for j, coef in enumerate(sol):
                        if coef is not None and coef != 0:
                            terms.append(f"{frac_to_latex(coef)} a_{{{j+1}}}")
                    if terms:
                        st.latex(r"b = " + " + ".join(terms))
                        st.write("Sustitución (vector reconstruido):")
                        st.latex(latex_column_vector(b_reconstruido))
                    else:
                        st.write("b = 0 (combinación trivial)")
                else:
                    st.success("✅ b pertenece al subespacio generado por las columnas de A (sistema consistente).")
                    st.subheader("⚡ Sistema con variables libres — solución parametrizada")
                    soluciones = ["0"] * (m-1)
                    for i in range(n):
                        piv_col = next((j for j in range(m-1) if R[i][j] != 0), None)
                        if piv_col is not None:
                            rhs = frac_str(R[i][-1])
                            expr = []
                            for j in free_vars:
                                if R[i][j] != 0:
                                    expr.append(f"-({frac_str(R[i][j])})*x{j+1}")
                            soluciones[piv_col] = rhs + (" + " + " + ".join(expr) if expr else "")
                    for j in free_vars:
                        soluciones[j] = f"x{j+1} (libre)"
                    for idx, val in enumerate(soluciones):
                        st.write(f"x{idx+1} = {val}")

            st.markdown('</div>', unsafe_allow_html=True)

# -------------------------
# SECCIÓN: AX = B
# -------------------------
elif seccion == "Ecuación matricial AX = B":
    st.header("🟰 Resolver AX = B (mejorado y verificación exacta)")

    A_text = st.text_area("Matriz A (filas separadas por saltos de línea)", value="1 2\n-2 5\n-5 6", height=140)
    B_text = st.text_input("Vector B (ej: 7 4 -3)", value="7 4 -3")

    if st.button("Ver resultado"):
        A_rows = []
        B_parts = []
        parse_error = False
        try:
            for row in A_text.splitlines():
                if row.strip():
                    parts = [p.strip() for p in row.replace(",", " ").split() if p.strip()]
                    fila = [Fraction(p) for p in parts]
                    A_rows.append(fila)
            B_parts = [Fraction(p) for p in B_text.replace(",", " ").split() if p.strip()]
        except Exception:
            st.error("⚠️ Error al parsear A o B. Usa enteros o fracciones (ej. -3, 4/5).")
            parse_error = True

        if parse_error or not A_rows:
            st.stop()

        if any(len(row) != len(A_rows[0]) for row in A_rows):
            st.error("⚠️ Todas las filas de A deben tener la misma cantidad de columnas.")
            st.stop()

        n_filas_A = len(A_rows)
        n_cols_A = len(A_rows[0])
        n_elem_B = len(B_parts)

        if n_elem_B < n_filas_A:
            B_parts = B_parts + [Fraction(0)] * (n_filas_A - n_elem_B)
        elif n_elem_B > n_filas_A:
            for _ in range(n_elem_B - n_filas_A):
                A_rows.append([Fraction(0)] * n_cols_A)
            n_filas_A = len(A_rows)

        M_aug = [ [Fraction(x) for x in row] + [Fraction(B_parts[i])] for i, row in enumerate(A_rows)]

        st.markdown('<div class="resultado">', unsafe_allow_html=True)
        st.latex(r"A = " + latex_matrix(A_rows))
        st.latex(r"B = " + latex_column_vector(B_parts))

        # Mostrar sistema de ecuaciones (como en tu foto)
        st.subheader("📑 Sistema de ecuaciones lineales (forma expandida)")
        eqs = pretty_equations(A_rows, B_parts, var_symbol="x")
        for e in eqs:
            st.latex(e)

        st.subheader("📋 Pasos de resolución (Gauss-Jordan)")
        pasos = gauss_jordan_steps(M_aug)
        if pasos:
            for i, (op, mat) in enumerate(pasos, start=1):
                st.markdown(f"**Paso {i}: {op}**")
                st.latex(latex_matrix_highlight(mat))
            R = pasos[-1][1]
        else:
            R = [row.copy() for row in M_aug]
            st.write("No se realizaron pasos (la matriz pudo estar ya en RREF).")
            st.latex(latex_matrix(R))

        st.subheader("✅ Matriz final (RREF):")
        st.latex(latex_matrix(R))

        # Interpretación
        n = len(R)
        m = len(R[0])
        pivot_cols = []
        pivot_row_for_col = {}
        for i in range(n):
            piv = None
            for j in range(m-1):
                if R[i][j] != 0:
                    piv = j
                    break
            if piv is not None:
                pivot_cols.append(piv)
                pivot_row_for_col[piv] = i
        pivot_cols = sorted(set(pivot_cols))
        free_vars = [j for j in range(m-1) if j not in pivot_cols]

        if not free_vars:
            sol = [None] * (m-1)
            for j in pivot_cols:
                row_idx = pivot_row_for_col[j]
                sol[j] = R[row_idx][-1]
            st.subheader("📌 Solución única encontrada:")
            for idx, val in enumerate(sol):
                if val is None:
                    st.write(f"x{idx+1} = libre")
                else:
                    st.write(f"x{idx+1} = {frac_str(val)}")

            columnas_A = []
            for col_idx in range(n_cols_A):
                col = [A_rows[row_idx][col_idx] for row_idx in range(n_filas_A)]
                columnas_A.append(col)

            b_reconstruido = [Fraction(0)] * n_filas_A
            for j, coef in enumerate(sol):
                if coef is not None:
                    for i in range(n_filas_A):
                        b_reconstruido[i] += coef * columnas_A[j][i]

            st.markdown("**Combinación lineal encontrada:**")
            terms = []
            for j, coef in enumerate(sol):
                if coef is not None and coef != 0:
                    terms.append(f"{frac_to_latex(coef)} a_{{{j+1}}}")
            if terms:
                st.latex(r"b = " + " + ".join(terms))
            else:
                st.write("b = 0 (combinación trivial)")

            st.subheader("🔍 Verificación de igualdad A·x = b")
            st.write("b (original):")
            st.latex(latex_column_vector(B_parts))
            st.write("b (reconstruido con A y solución x):")
            st.latex(latex_column_vector(b_reconstruido))

            if b_reconstruido == [Fraction(x) for x in B_parts]:
                st.success("✅ Verificación correcta: A·x coincide exactamente con b.")
            else:
                st.error("❌ Verificación fallida: A·x NO coincide con b. Diferencias:")
                dif = [b_reconstruido[i] - Fraction(B_parts[i]) for i in range(n_filas_A)]
                for i, d in enumerate(dif, start=1):
                    st.write(f"Componente {i}: diferencia = {frac_str(d)}")
        else:
            st.subheader("⚡ Sistema con variables libres (solución parametrizada)")
            soluciones = ["0"] * (m-1)
            for i in range(n):
                piv_col = next((j for j in range(m-1) if R[i][j] != 0), None)
                if piv_col is not None:
                    rhs = frac_str(R[i][-1])
                    expr = []
                    for j in free_vars:
                        if R[i][j] != 0:
                            expr.append(f"-({frac_str(R[i][j])})*x{j+1}")
                    soluciones[piv_col] = rhs + (" + " + " + ".join(expr) if expr else "")
            for j in free_vars:
                soluciones[j] = f"x{j+1} (libre)"
            for idx, val in enumerate(soluciones):
                st.write(f"x{idx+1} = {val}")

            st.info("Como el sistema tiene variables libres, b pertenece al subespacio generado por las columnas de A si y sólo si el sistema es consistente (no hay filas del tipo [0 ... 0 | c] con c != 0).")

        st.markdown("### Columnas de A (vectores generadores)")
        col_vectors_latex = [latex_column_vector([A_rows[r][c] for r in range(n_filas_A)]) for c in range(n_cols_A)]
        for i, col_ltx in enumerate(col_vectors_latex, start=1):
            st.latex(r"a_" + str(i) + " = " + col_ltx)

        st.markdown('</div>', unsafe_allow_html=True)

# -------------------------
# SECCIÓN: GAUSS / Gauss-Jordan
# -------------------------
elif seccion == "Gauss / Gauss-Jordan (RREF)":
    st.header("⚙️ Gauss / Gauss-Jordan (entradas numéricas obligatorias)")
    n_filas = st.number_input("Número de Filas (matriz aumentada):", min_value=1, max_value=10, value=3, step=1, key="g_nfilas")
    n_columnas = st.number_input("Número de Columnas (incluye columna independiente):", min_value=1, max_value=12, value=4, step=1, key="g_ncols")

    st.write("Rellena la matriz aumentada. Deja en blanco para 0. **Sólo números/ fracciones** (p. ej. -1, 3/2).")
    matriz_inputs = []
    invalid_cells = []
    for i in range(n_filas):
        cols = st.columns(n_columnas)
        fila_vals = []
        for j in range(n_columnas):
            with cols[j]:
                val = st.text_input(f"({i+1},{j+1})", "0", key=f"gauss_{i}_{j}")
                sval = str(val).strip()
                if sval == "":
                    sval = "0"
                try:
                    frac = Fraction(sval)
                    fila_vals.append(frac)
                except Exception:
                    fila_vals.append(None)
                    invalid_cells.append((i+1, j+1, sval))
        matriz_inputs.append(fila_vals)

    if invalid_cells:
        st.error("Hay entradas inválidas en la matriz (deben ser números o fracciones). Ejemplo de celdas inválidas:")
        for (fi, cj, sval) in invalid_cells:
            st.write(f"Fila {fi}, Columna {cj}: '{sval}'")
    else:
        if st.button("Ver resultado"):
            st.markdown('<div class="resultado">', unsafe_allow_html=True)
            st.latex(latex_matrix(matriz_inputs))
            pasos = gauss_jordan_steps(matriz_inputs)
            if pasos:
                for k, (op, mat) in enumerate(pasos, start=1):
                    st.write(f"🔹 Paso {k}: {op}")
                    st.latex(latex_matrix_highlight(mat, []))
                st.subheader("✅ Matriz final (RREF):")
                st.latex(latex_matrix(pasos[-1][1]))
                resultado_texto, pivot_cols, free_vars = resolver_sistema(matriz_inputs)
                st.subheader("📋 Análisis de solución:")
                st.code(resultado_texto, language="text")
                st.markdown("### Explicación de variables")
                st.write("Variables básicas (columnas/pivote):", [f"x{p+1}" for p in pivot_cols] if pivot_cols else "Ninguna")
                st.write("Variables libres:", [f"x{f+1}" for f in free_vars] if free_vars else "Ninguna")
            else:
                st.write("No se realizaron pasos (la matriz ya estaba en RREF o no se encontró pivote).")
                st.latex(latex_matrix(matriz_inputs))
            st.markdown('</div>', unsafe_allow_html=True)

# -------------------------
# Footer
# -------------------------
st.markdown("---")
st.write("🔢 ")