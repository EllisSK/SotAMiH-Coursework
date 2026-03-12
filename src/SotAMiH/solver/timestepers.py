def calculate_timestep(cells, resolution):
    g = 9.81
    c = 1.0

    max_speed = 0

    for cell in cells:
        h = cell.Q_vector[0]
        hu = cell.Q_vector[1]

        if h > 1e-6:
            u = hu / h

            local_speed = float((abs(u) + ((g * h) ** (0.5))))

            if local_speed > max_speed:
                max_speed = local_speed

    if max_speed == 0:
        return 0.1
    
    dt = c * resolution / max_speed

    return dt

def first_order_time_marching(F_left, Q_vector, F_right, S_vector, timestep, resolution):
    return Q_vector - ((timestep / resolution) * (F_right - F_left)) + (timestep * S_vector)

def rk2_time_marching():
    pass

def rk4_time_marching():
    pass