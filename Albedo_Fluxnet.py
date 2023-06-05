"""Theoratical"""

MIN_ALB = 0.5
MAX_ALB = 0.85

COLD_SNOW_RELAX = 0.008
MELT_RELAX = 0.24

timerange = np.arange(0,15)
snowfall_event = [3, 8]
first_alb = float(0.75)


def surface_albedo_to_clear_sky(Albedo):
    # Chen et al 1983
    return(0.0587+0.730*Albedo)

def Albedo_Scheme(MIN_ALB, MAX_ALB, COLD_RELAX, WARM_RELAX, SNOWFALL, TEMPERATURE):



def THEORY_ALB_NORMAL(MIN_ALB, MAX_ALB, start_alb, RELAX, timerange):
    albedo = np.zeros(len(timerange))
    albedo[0] = start_alb
    for time in timerange[0:-1]:
        albedo[time+1] = albedo[time] - RELAX
        if albedo[time+1] < MIN_ALB:
            albedo[time+1] = MIN_ALB
    return(albedo)
def THEORY_ALB_SNOWFALL(MIN_ALB, MAX_ALB, start_alb, RELAX, timerange, snowfall_time):
    albedo = np.zeros(len(timerange))
    albedo[0] = start_alb
    for time in timerange[0:-1]:
        albedo[time+1] = albedo[time] - RELAX
        if albedo[time+1] < MIN_ALB:
            albedo[time+1] = MIN_ALB

        if time in snowfall_time:
            albedo[time+1] = MAX_ALB
    return(albedo)
def THEORY_ALB_MELT(MIN_ALB, MAX_ALB, start_alb, RELAX, time, snowfall_time):
    albedo = np.zeros(len(timerange))
    albedo[0] = start_alb
    for time in timerange[0:-1]:
        albedo[time+1] = (albedo[time] - MIN_ALB) * np.exp(-RELAX) +MIN_ALB
        if albedo[time+1] < MIN_ALB:
            albedo[time+1] = MIN_ALB
        if time in snowfall_time:
            albedo[time+1] = MAX_ALB
    return(albedo)