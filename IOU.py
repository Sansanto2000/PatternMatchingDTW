def IoU(obs_min:float, obs_max:float, teo_min:float, teo_max:float) -> float:
    """Funcion para realizar el calculo de [Interseccion / Union] dados dos rangos de valores

    Args:
        obs_min (float): Valor minimo del rango 1
        obs_max (float): Valor maximo del rango 1
        teo_min (float): Valor minimo del rango 2
        teo_max (float): Valor maximo del rango 2

    Returns:
        float: Interseccion/Union
    """

    if (obs_min > teo_max or teo_min > obs_max):
        interseccion = 0
        # Obs:        |----| 
        # Teo: |----|
        #      or
        # Obs: |----| 
        # Teo:        |----|
        
    elif (teo_min <= obs_min and teo_max >= obs_max):
        interseccion = obs_max - obs_min 
        # Obs:   |----| 
        # Teo: |--------|
        
    elif (teo_min >= obs_min and teo_max <= obs_max):
        interseccion = teo_max - teo_min
        # Obs: |--------| 
        # Teo:   |----|
        
    elif (obs_min <= teo_min and teo_min <= obs_max and obs_max <= teo_max):
        interseccion = obs_max - teo_min
        # Obs: |----| 
        # Teo:    |----|
        
    elif (teo_min <= obs_min and obs_min <= teo_max and teo_max <= obs_max):
        interseccion = teo_max - obs_min
        # Obs:    |----| 
        # Teo: |----|

    union = teo_max - teo_min + obs_max - obs_min - interseccion
    iou = union / interseccion
    return iou