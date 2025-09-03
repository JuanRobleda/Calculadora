from flask import Flask, render_template, request, jsonify
import webbrowser
import threading

app = Flask(__name__)

# ====== FUNCIONES DE MATRICES ======
def format_matrix(M):
    return "\n".join(["\t".join(f"{num:.2f}" for num in fila) for fila in M])

def determinante_rec(M):
    if len(M) == 1:
        return M[0][0]
    if len(M) == 2:
        return M[0][0]*M[1][1] - M[0][1]*M[1][0]
    det = 0
    for c in range(len(M)):
        sub = [fila[:c] + fila[c+1:] for fila in (M[1:])]
        det += ((-1)**c) * M[0][c] * determinante_rec(sub)
    return det

def gauss_jordan_proceso(M):
    n = len(M)
    m = len(M[0])
    A = [fila[:] for fila in M]
    pasos = ["Matriz inicial:\n" + format_matrix(A)]

    for i in range(n):
        if A[i][i] == 0:
            for k in range(i+1, n):
                if A[k][i] != 0:
                    A[i], A[k] = A[k], A[i]
                    pasos.append(f"Intercambio fila {i+1} con fila {k+1}:\n" + format_matrix(A))
                    break
        pivote = A[i][i]
        if pivote == 0:
            raise Exception("No se puede continuar, pivote cero.")

        A[i] = [x / pivote for x in A[i]]
        pasos.append(f"Normalizar fila {i+1}:\n" + format_matrix(A))

        for j in range(n):
            if i != j:
                factor = A[j][i]
                A[j] = [A[j][k] - factor*A[i][k] for k in range(m)]
        pasos.append(f"Eliminación en columna {i+1}:\n" + format_matrix(A))

    return pasos, A

# ====== RUTAS FLASK ======
@app.route("/")
def index():
    return render_template("index.html")

@app.route("/mostrar", methods=["POST"])
def mostrar():
    try:
        data = request.json
        matriz = [[float(x) for x in fila] for fila in data["matriz"]]
        return jsonify({"ok": True, "resultado": format_matrix(matriz)})
    except Exception as e:
        return jsonify({"ok": False, "error": str(e)})

@app.route("/gauss", methods=["POST"])
def gauss():
    try:
        data = request.json
        matriz = [[float(x) for x in fila] for fila in data["matriz"]]
        pasos, reducida = gauss_jordan_proceso(matriz)
        texto = "\n\n".join(pasos)
        texto += "\n\n✅ Matriz final:\n" + format_matrix(reducida)
        return jsonify({"ok": True, "resultado": texto})
    except Exception as e:
        return jsonify({"ok": False, "error": str(e)})

# ====== FUNCION PARA ABRIR NAVEGADOR ======
def abrir_navegador():
    webbrowser.open_new("http://127.0.0.1:5000/")

if __name__ == "__main__":
    # Ejecutar Flask en un hilo separado para abrir navegador al mismo tiempo
    threading.Timer(1.5, abrir_navegador).start()
    app.run(debug=True)