def get_simulation_params(pde_model='volodin_2025'):
    if pde_model == 'volodin_2025':
        return {
            't_range': 500,
            'dt': 1e-4,
            'Ca': 0.00072,
            'G0': 0.0000098,
            'phi': 1e-6,
            'Pr': 6.99,
            'Bi': 1,
            'Ma': 0.167,
            'omega': 1.5,
            'b': 2,
            'bc': {'x': 'periodic'},
        }
    elif pde_model == 'shklyaev_2008_viii_b':
        return {
            't_range': 500,
            'V': 0,
            'k': 2.2,
            'G0': 3.33 * 10**-4,
            'omega': 0.2,
            'bc': {'x': 'periodic'},
        }
    else:
        raise ValueError(f'Unknown model_name: {pde_model}')
