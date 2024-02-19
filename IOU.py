def IoU(obs_min:float, obs_max:float, real_min:float, real_max:float) -> float:
    """Funcion para realizar el calculo de [Interseccion / Union] dados dos rangos de valores

    Args:
        obs_min (float): Valor minimo del rango 1
        obs_max (float): Valor maximo del rango 1
        real_min (float): Valor minimo del rango 2
        real_max (float): Valor maximo del rango 2

    Returns:
        float: Interseccion/Union
    """

    if (obs_min > real_max or real_min > obs_max):
        interseccion = 0
        # Obs:        |----| 
        # Rea: |----|
        #      or
        # Obs: |----| 
        # Rea:        |----|
        
    elif (real_min <= obs_min and real_max >= obs_max):
        interseccion = obs_max - obs_min 
        # Obs:   |----| 
        # Rea: |--------|
        
    elif (real_min >= obs_min and real_max <= obs_max):
        interseccion = real_max - real_min
        # Obs: |--------| 
        # Rea:   |----|
        
    elif (obs_min <= real_min and real_min <= obs_max and obs_max <= real_max):
        interseccion = obs_max - real_min
        # Obs: |----| 
        # Rea:    |----|
        
    elif (real_min <= obs_min and obs_min <= real_max and real_max <= obs_max):
        interseccion = real_max - obs_min
        # Obs:    |----| 
        # Rea: |----|

    union = real_max - real_min + obs_max - obs_min - interseccion
    iou =  interseccion / union
    return iou