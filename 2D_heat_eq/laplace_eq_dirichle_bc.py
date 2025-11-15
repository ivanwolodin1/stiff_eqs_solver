from pde import CartesianGrid, solve_laplace_equation

grid = CartesianGrid([[0, 5], [0, 5]], 100)
bcs = {
    'x+': {'value': '0'},
    'x-': {'value': '0'},
    'y+': {'value': 'sin(x)'},
    'y-': {'value': '0'},
}

# what solver is used here?
res = solve_laplace_equation(grid, bcs)
res.plot()
