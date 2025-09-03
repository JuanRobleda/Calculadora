import tkinter as tk
from tkinter import messagebox

class MatrizApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Operaciones con Matrices")

        tk.Label(root, text="Número de Filas:").pack()
        self.filas_entry = tk.Entry(root)
        self.filas_entry.pack()

        tk.Label(root, text="Número de Columnas:").pack()
        self.columnas_entry = tk.Entry(root)
        self.columnas_entry.pack()

        
        self.crear_btn = tk.Button(root, text="Crear Matriz", command=self.crear_matriz)
        self.crear_btn.pack()

        self.entries = []

    def crear_matriz(self):
    
        for fila in self.entries:
            for e in fila:
                e.destroy()
        self.entries.clear()

        try:
            filas = int(self.filas_entry.get())
            columnas = int(self.columnas_entry.get())
        except:
            messagebox.showerror("Error", "Introduce enteros válidos")
            return

        self.matriz_frame = tk.Frame(self.root)
        self.matriz_frame.pack()

        for i in range(filas):
            fila_entries = []
            for j in range(columnas):
                e = tk.Entry(self.matriz_frame, width=5, justify="center")
                e.grid(row=i, column=j, padx=2, pady=2)
                fila_entries.append(e)
            self.entries.append(fila_entries)

        
        tk.Button(self.root, text="Mostrar Matriz", command=self.mostrar_matriz).pack()
        tk.Button(self.root, text="Gauss-Jordan", command=self.gauss_jordan).pack()

    def leer_matriz(self):
        try:
            return [[float(self.entries[i][j].get()) for j in range(len(self.entries[0]))] 
                    for i in range(len(self.entries))]
        except:
            messagebox.showerror("Error", "Introduce solo números válidos")
            return None

    def mostrar_matriz(self):
        A = self.leer_matriz()
        if A:
            messagebox.showinfo("Matriz", self.format_matrix(A))


    def gauss_jordan(self):
        A = self.leer_matriz()
        if A:
            try:
                pasos, reducida = self.gauss_jordan_proceso(A)
                texto = "\n\n".join(pasos)
                texto += "\n\n✅ Matriz final:\n" + self.format_matrix(reducida)
                messagebox.showinfo("Gauss-Jordan - Procedimiento", texto)
            except Exception as e:
                messagebox.showerror("Error", str(e))

    def format_matrix(self, M):
        return "\n".join(["\t".join(f"{num:.2f}" for num in fila) for fila in M])

    def determinante_rec(self, M):
        if len(M) == 1:
            return M[0][0]
        if len(M) == 2:
            return M[0][0]*M[1][1] - M[0][1]*M[1][0]
        det = 0
        for c in range(len(M)):
            sub = [fila[:c] + fila[c+1:] for fila in (M[1:])]
            det += ((-1)**c) * M[0][c] * self.determinante_rec(sub)
        return det

    def gauss_jordan_proceso(self, M):
        n = len(M)
        m = len(M[0])
        A = [fila[:] for fila in M]
        pasos = ["Matriz inicial:\n" + self.format_matrix(A)]

        for i in range(n):
            
            if A[i][i] == 0:
                for k in range(i+1, n):
                    if A[k][i] != 0:
                        A[i], A[k] = A[k], A[i]
                        pasos.append(f"Intercambio fila {i+1} con fila {k+1}:\n" + self.format_matrix(A))
                        break
            pivote = A[i][i]
            if pivote == 0:
                raise Exception("No se puede continuar, pivote cero.")

            A[i] = [x / pivote for x in A[i]]
            pasos.append(f"Normalizar fila {i+1}:\n" + self.format_matrix(A))

            
            for j in range(n):
                if i != j:
                    factor = A[j][i]
                    A[j] = [A[j][k] - factor*A[i][k] for k in range(m)]
            pasos.append(f"Eliminación en columna {i+1}:\n" + self.format_matrix(A))

        return pasos, A

    def gauss_jordan_inverse(self, M):
        n = len(M)
        A = [fila[:] + [1 if i == j else 0 for j in range(n)] for i, fila in enumerate(M)]
        _, R = self.gauss_jordan_proceso(A)
        return [fila[n:] for fila in R]


if __name__ == "__main__":
    root = tk.Tk()
    app = MatrizApp(root)
    root.mainloop()