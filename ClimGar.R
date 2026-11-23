## ClimGar — Clasificación climática de García (sistema mexicano)
##   bd  : data.frame con columnas Station, Month (1-12), Tc (temp °C), Pmm (prec mm)
##   y   : hemisferio forzado "N" o "S"; si es NA se infiere automáticamente
##   - Invierno   = 3 meses consecutivos más fríos   (ventana circular sobre meses)
##   - Verano     = 6 meses consecutivos más cálidos con desfase +1 mes adelante

ClimGar <- function(bd, y = NA) {
  
  require(dplyr)
  require(stringr)

  # ── 0. Estaciones con los 12 meses completos ──────────────────────────────
  estaciones_ok <- bd |>
    drop_na() |>
    group_by(Station) |>
    summarise(nn = n_distinct(Month), .groups = "drop") |>
    dplyr::filter(nn == 12) |>
    pull(Station)
  
  if (length(estaciones_ok) == 0) {
    return(data.frame(
      Station=NA,Pmm=NA,Pmin=NA,Pmax=NA,Tmean=NA,Tmin=NA,Tmax=NA,
      Pver=NA,Regi=NA,P_T=NA,Osc=NA,PPInv=NA,Hum=NA,Hemi=NA,
      Tclim1=NA,Tclim2=NA,Tclim3=NA,Tclim4=NA,Tclim5=NA,Ampli=NA,Clim=NA
    )[-1, ])
  }
  
  bd_ok <- bd |> dplyr::filter(Station %in% estaciones_ok) |> drop_na()
  
  # ── 1. Estadísticos básicos por estación (todo vectorizado) ───────────────
  base <- bd_ok |>
    group_by(Station) |>
    summarise(
      Ptot  = sum(Pmm),
      Pmin  = min(Pmm),
      Pmax  = max(Pmm),
      Tmean = round(mean(Tc), 2),
      Tmin  = min(Tc),
      Tmax  = max(Tc),
      md    = sum(Tc > 10),          # meses con T>10°C
      .groups = "drop"
    ) |>
    mutate(
      P_T = round(Ptot / Tmean, 2),
      Osc = Tmax - Tmin
    )
  
  # ── 2. Hemisferio ─────────────────────────────────────────────────────────
  # Para inferirlo: comparamos máxima Tc en meses cálidos boreales (abr-jul)
  # vs máxima Tc en meses cálidos australes (oct-ene)
  hemi_infer <- bd_ok |>
    group_by(Station) |>
    summarise(
      Tc_boreal  = max(Tc[Month %in% 4:7]),
      Tc_austral = max(Tc[Month %in% c(10:12, 1)]),
      Osc_loc    = max(Tc) - min(Tc),
      .groups = "drop"
    ) |>
    mutate(
      Hemi = case_when(
        !is.na(!!y) ~ y,                          # forzado por el usuario
        Osc_loc < 1  ~ "E",                       # ecuatorial, sin estacionalidad
        Tc_boreal >= Tc_austral ~ "N",
        TRUE ~ "S"
      )
    ) |>
    dplyr::select(Station, Hemi)
  
  # Si y está forzado, sobreescribir todo
  if (!is.na(y)) {
    hemi_infer <- hemi_infer |> mutate(Hemi = y)
  }
  
  base <- base |> left_join(hemi_infer, by = "Station")
  
  # ── 3. Detección circular de ventanas de invierno y verano ────────────────
  #
  # Invierno = inicio del bloque de 3 meses consecutivos con menor suma de Tc
  # Verano   = inicio del bloque de 6 meses consecutivos con mayor suma de Tc,
  #            desplazado +1 mes adelante (desfase térmico océano-atmósfera)
  #
  # Trabajamos con índice circular: mes duplicado 1-12-1-12 para ventanas
  
  temp_mensual <- bd_ok |>
    group_by(Station, Month) |>
    summarise(Tc_m = mean(Tc), Pmm_m = mean(Pmm), .groups = "drop")
  
  # Función auxiliar: suma de ventana circular de longitud `w` empezando en mes `s`
  suma_circ <- function(vec, s, w) {
    idx <- ((s - 1 + 0:(w-1)) %% 12) + 1
    sum(vec[idx])
  }
  
  # Para cada estación, encontrar mes de inicio de invierno y verano
  season_windows <- temp_mensual |>
    group_by(Station) |>
    summarise(
      # Vector ordenado por mes 1-12
      Tc_vec  = list(Tc_m[order(Month)]),
      Pmm_vec = list(Pmm_m[order(Month)]),
      .groups = "drop"
    ) |>
    rowwise() |>
    mutate(
      # Suma de Tc para cada ventana de 3 meses consecutivos
      sumas_inv = list(sapply(1:12, \(s) suma_circ(Tc_vec, s, 3))),
      # Mes de inicio del invierno (bloque más frío)
      inv_start = which.min(sumas_inv),
      # Meses de invierno (3 consecutivos a partir de inv_start)
      inv_months = list(((inv_start - 1 + 0:2) %% 12) + 1),
      
      # Suma de Tc para cada ventana de 6 meses consecutivos
      sumas_ver = list(sapply(1:12, \(s) suma_circ(Tc_vec, s, 6))),
      # Mes de inicio del verano = bloque más cálido, +1 mes de desfase
      ver_start_raw = which.max(sumas_ver),
      ver_start = (ver_start_raw %% 12) + 1,   # desfase +1 mes
      # Meses de verano (6 consecutivos)
      ver_months = list(((ver_start - 1 + 0:5) %% 12) + 1),
      
      # Meses que NO son ni verano ni invierno → transición
      trans_months = list(setdiff(1:12, union(inv_months, ver_months)))
    ) |>
    ungroup()
  
  # ── 4. Calcular Pver, PPInv, Pmav, Pmah, Pmiv por estación ───────────────
  season_stats <- season_windows |>
    left_join(
      temp_mensual |>
        group_by(Station) |>
        summarise(Pmm_list = list(setNames(Pmm_m[order(Month)], 1:12)), .groups="drop"),
      by = "Station"
    ) |>
    rowwise() |>
    mutate(
      Pver  = sum(Pmm_list[ver_months]),
      Pinv  = sum(Pmm_list[inv_months]),
      Pmav  = max(Pmm_list[ver_months]),     # max mensual en verano
      Pmiv  = min(Pmm_list[ver_months]),     # min mensual en verano
      Pmah  = max(Pmm_list[inv_months]),     # max mensual en invierno
      PPInv = round(Pinv / sum(Pmm_list) * 100, 2),
      # ¿Hay canícula? (mínimo local en pleno verano, meses 7-8 para NH)
      # Se evalúa después de conocer Hemi; aquí precalculamos los valores crudos
      Pmm7  = Pmm_list[[7]],
      Pmm6  = Pmm_list[[6]],
      Pmm8  = Pmm_list[[8]],
      Pmm1  = Pmm_list[[1]],
      Pmm2  = Pmm_list[[2]],
      Pmm12 = Pmm_list[[12]],
      Pmm3  = Pmm_list[[3]]
    ) |>
    ungroup() |>
    select(Station, Pver, Pinv, Pmav, Pmiv, Pmah, PPInv,
           Pmm7, Pmm6, Pmm8, Pmm1, Pmm2, Pmm12, Pmm3)
  
  base <- base |>
    left_join(season_stats, by = "Station")
  
  # ── 5. Régimen de lluvias ──────────────────────────────────────────────────
  # ver  = verano lluvioso: Pmav > 10×Pmin
  # inv  = invierno lluvioso: Pmah > 3×Pmiv AND Pmax == Pmah
  # int  = intermedio
  
  base <- base |>
    mutate(
      Regi = case_when(
        Pmav > 10 * Pmin                    ~ "ver",
        Pmah > 3  * Pmiv & Pmax == Pmah    ~ "inv",
        TRUE                                 ~ "int"
      )
    )
  
  # ── 6. Amplitud térmica (Tclim5 / Ampli) ─────────────────────────────────
  base <- base |>
    mutate(
      Ampli = case_when(
        Osc > 14 ~ "(e')",
        Osc > 7  ~ "(e)",
        Osc > 5  ~ "(i')",
        TRUE     ~ "i"
      )
    )
  
  # ── 7. Tipo climático principal — Tclim1 ──────────────────────────────────
  base <- base |>
    mutate(
      # Primero clasificar áridos/semiáridos (B), luego los demás
      Tclim1 = case_when(
        # B: áridos — umbral depende del régimen y PPInv
        Regi == "ver" & PPInv <  10.2 & Ptot/10 < Tmean + 14  ~ "BW",
        Regi == "ver" & PPInv <  10.2 & Ptot/10 < 2*Tmean+28 & P_T >= 22.9 ~ "BS1",
        Regi == "ver" & PPInv <  10.2 & Ptot/10 < 2*Tmean+28 & P_T <  22.9 ~ "BS0",
        Regi == "ver" & PPInv >= 10.2 & Ptot/10 < Tmean + 10.5 ~ "BW",
        Regi == "ver" & PPInv >= 10.2 & Ptot/10 < 2*Tmean+21 & P_T >= 22.9 ~ "BS1",
        Regi == "ver" & PPInv >= 10.2 & Ptot/10 < 2*Tmean+21 & P_T <  22.9 ~ "BS0",
        Regi == "int" & Ptot/10 < Tmean + 7   ~ "BW",
        Regi == "int" & Ptot/10 < 2*Tmean+14 & P_T >= 22.9 ~ "BS1",
        Regi == "int" & Ptot/10 < 2*Tmean+14 & P_T <  22.9 ~ "BS0",
        Regi == "inv" & Ptot/10 < Tmean       ~ "BW",
        Regi == "inv" & Ptot/10 < 2*Tmean    & P_T >= 22.9 ~ "BS1",
        Regi == "inv" & Ptot/10 < 2*Tmean    & P_T <  22.9 ~ "BS0",
        # E: polar
        Tmax < 10                             ~ "E",
        # D: continental frío
        Tmax >= 10 & Tmin < -3               ~ "D",
        # A: cálido húmedo (evaluar antes que C — C no tiene cota superior
        # de Tmean y de lo contrario nunca se alcanzan estas ramas)
        Tmean > 22 & Tmin >  18             ~ "A",
        Tmean > 22 & Tmin <= 18             ~ "A(C)",
        # Transición (A)C
        Tmean > 18 & Tmean <= 22             ~ "(A)C",
        # C: templado
        Tmean > 5  & Tmin >= -3              ~ "C",
        TRUE ~ NA_character_
      )
    )
  
  # ── 8. Segunda letra — Tclim2 (régimen hídrico) ───────────────────────────
  base <- base |>
    mutate(
      # Humedad de fondo
      graph = case_when(
        Tclim1 %in% c("A","A(C)") & (Pmin > 60 | Pmin*25 < Ptot-1000) ~ "hum",
        Tclim1 %in% c("A","A(C)")                                        ~ "subhum",
        Tclim1 %in% c("C","D","(A)C","E") & (Pmin > 40 | Pmin*32.5 < Ptot-500) ~ "hum",
        Tclim1 %in% c("C","D","(A)C","E")                               ~ "subhum",
        TRUE ~ NA_character_
      ),
      
      Tclim2 = case_when(
        # ── RÉGIMEN VERANO ──
        Regi == "ver" & Tclim1 %in% c("A","A(C)") & graph == "hum" & PPInv < 5    ~ "m(w)",
        Regi == "ver" & Tclim1 %in% c("A","A(C)") & graph == "hum" & PPInv < 10.2 ~ "m",
        Regi == "ver" & Tclim1 %in% c("A","A(C)") & graph == "hum"                ~ "m(f)",
        
        Regi == "ver" & Tclim1 %in% c("C","D","(A)C") & graph == "hum" & PPInv < 5    ~ "(m)(w)",
        Regi == "ver" & Tclim1 %in% c("C","D","(A)C") & graph == "hum" & PPInv < 10.2 ~ "(m)",
        Regi == "ver" & Tclim1 %in% c("C","D","(A)C") & graph == "hum"                ~ "(m)(f)",
        
        Regi == "ver" & Tclim1 %in% c("C","D","(A)C","E") & graph == "subhum" & P_T > 55  & PPInv < 5    ~ "(w2)(w)",
        Regi == "ver" & Tclim1 %in% c("C","D","(A)C","E") & graph == "subhum" & P_T > 55  & PPInv < 10.2 ~ "(w2)",
        Regi == "ver" & Tclim1 %in% c("C","D","(A)C","E") & graph == "subhum" & P_T > 55                 ~ "(w2)(x')",
        
        Regi == "ver" & Tclim1 %in% c("A","A(C)") & graph == "subhum" & P_T >= 55.3 & PPInv < 5    ~ "w2(w)",
        Regi == "ver" & Tclim1 %in% c("A","A(C)") & graph == "subhum" & P_T >= 55.3 & PPInv < 10.2 ~ "w2",
        Regi == "ver" & Tclim1 %in% c("A","A(C)") & graph == "subhum" & P_T >= 55.3               ~ "w2(x')",
        
        Regi == "ver" & Tclim1 %in% c("C","D","(A)C","E") & graph == "subhum" & P_T < 55 & P_T >= 43.2 & Pmin < (100-Ptot/25) & PPInv < 5    ~ "(w1)(w)",
        Regi == "ver" & Tclim1 %in% c("C","D","(A)C","E") & graph == "subhum" & P_T < 55 & P_T >= 43.2 & Pmin < (100-Ptot/25) & PPInv < 10.2 ~ "(w1)",
        Regi == "ver" & Tclim1 %in% c("C","D","(A)C","E") & graph == "subhum" & P_T < 55 & P_T >= 43.2 & Pmin < (100-Ptot/25)               ~ "(w1)(x')",
        
        Regi == "ver" & Tclim1 %in% c("A","A(C)") & graph == "subhum" & P_T < 55.3 & P_T > 43.2 & PPInv < 5    ~ "w1(w)",
        Regi == "ver" & Tclim1 %in% c("A","A(C)") & graph == "subhum" & P_T < 55.3 & P_T > 43.2 & PPInv < 10.2 ~ "w1",
        Regi == "ver" & Tclim1 %in% c("A","A(C)") & graph == "subhum" & P_T < 55.3 & P_T > 43.2               ~ "w1(x')",
        
        Regi == "ver" & Tclim1 %in% c("C","D","(A)C","E") & graph == "subhum" & P_T < 43.2 & PPInv < 5    ~ "(w0)(w)",
        Regi == "ver" & Tclim1 %in% c("C","D","(A)C","E") & graph == "subhum" & P_T < 43.2 & PPInv < 10.2 ~ "(w0)",
        Regi == "ver" & Tclim1 %in% c("C","D","(A)C","E") & graph == "subhum" & P_T < 43.2               ~ "(w0)(x')",
        
        Regi == "ver" & Tclim1 %in% c("A","A(C)") & graph == "subhum" & P_T <= 43.2 & PPInv < 5    ~ "w0(w)",
        Regi == "ver" & Tclim1 %in% c("A","A(C)") & graph == "subhum" & P_T <= 43.2 & PPInv < 10.2 ~ "w0",
        Regi == "ver" & Tclim1 %in% c("A","A(C)") & graph == "subhum" & P_T <= 43.2               ~ "w0(x')",
        
        Regi == "ver" & Tclim1 %in% c("BS0","BS1","BW") & PPInv < 5    ~ "w(w)",
        Regi == "ver" & Tclim1 %in% c("BS0","BS1","BW") & PPInv < 10.2 ~ "w",
        Regi == "ver" & Tclim1 %in% c("BS0","BS1","BW")                ~ "w(x')",
        
        # ── RÉGIMEN INTERMEDIO ──
        Regi == "int" & Tclim1 %in% c("A","A(C)") & graph == "hum" & PPInv < 18 ~ "f(m)",
        Regi == "int" & Tclim1 %in% c("A","A(C)") & graph == "hum"               ~ "f",
        
        Regi == "int" & Tclim1 %in% c("C","D","(A)C","E") & graph == "hum" & PPInv < 18 ~ "(fm)",
        Regi == "int" & Tclim1 %in% c("C","D","(A)C","E") & graph == "hum"               ~ "(f)",
        
        Regi == "int" & Tclim1 %in% c("C","D","(A)C","E") & graph == "subhum" & P_T > 55    & PPInv < 18 ~ "(x')(w2)",
        Regi == "int" & Tclim1 %in% c("C","D","(A)C","E") & graph == "subhum" & P_T > 55                 ~ "x'",
        
        Regi == "int" & Tclim1 %in% c("A","A(C)") & graph == "subhum" & P_T >= 55.3 & PPInv < 18 ~ "x'(w2)",
        Regi == "int" & Tclim1 %in% c("A","A(C)") & graph == "subhum" & P_T >= 55.3               ~ "x'",
        
        Regi == "int" & Tclim1 %in% c("C","D","(A)C","E") & graph == "subhum" & P_T < 55 & P_T >= 43.2 & Pmin < (100-Ptot/25) & PPInv < 18 ~ "(x')(w1)",
        Regi == "int" & Tclim1 %in% c("C","D","(A)C","E") & graph == "subhum" & P_T < 55 & P_T >= 43.2 & Pmin < (100-Ptot/25)               ~ "x'",
        
        Regi == "int" & Tclim1 %in% c("A","A(C)") & graph == "subhum" & P_T < 55.3 & P_T >= 43.2 & PPInv < 18 ~ "x'(w1)",
        Regi == "int" & Tclim1 %in% c("A","A(C)") & graph == "subhum" & P_T < 55.3 & P_T >= 43.2               ~ "x'",
        
        Regi == "int" & Tclim1 %in% c("C","D","(A)C","E") & graph == "subhum" & P_T < 43.2 & PPInv < 18 ~ "(x')(w0)",
        Regi == "int" & Tclim1 %in% c("C","D","(A)C","E") & graph == "subhum" & P_T < 43.2               ~ "x'",
        
        Regi == "int" & Tclim1 %in% c("A","A(C)") & graph == "subhum" & P_T < 43.2 & PPInv < 18 ~ "x'(w0)",
        Regi == "int" & Tclim1 %in% c("A","A(C)") & graph == "subhum" & P_T < 43.2               ~ "x'",
        
        Regi == "int" & Tclim1 %in% c("BS0","BS1","BW") & PPInv < 18 ~ "x'(w)",
        Regi == "int" & Tclim1 %in% c("BS0","BS1","BW")               ~ "x'",
        
        # ── RÉGIMEN INVIERNO ──
        Regi == "inv" & (Pmin > 60 | Pmin > (100 - Ptot/25)) & PPInv < 36 ~ "f(s)",
        Regi == "inv" & (Pmin > 60 | Pmin > (100 - Ptot/25))               ~ "s(f)",
        Regi == "inv" & PPInv < 36 ~ "s(x')",
        Regi == "inv"              ~ "s",
        
        TRUE ~ NA_character_
      )
    )
  
  # ── 9. Tercera letra — Tclim3 (temperatura del mes más cálido) ────────────
  base <- base |>
    mutate(
      Tclim3 = case_when(
        # B: subcontinental h/k
        Tclim1 %in% c("BW","BS0","BS1") & Tmean > 22 & Tmin > 18  ~ "h'",
        Tclim1 %in% c("BW","BS0","BS1") & Tmean > 22               ~ "(h')h",
        Tclim1 %in% c("BW","BS0","BS1") & Tmean > 18               ~ "h'(h)",
        Tclim1 %in% c("BW","BS0","BS1") & Tmin > 18 & Tmax > 22   ~ "h",
        Tclim1 %in% c("BW","BS0","BS1") & Tmean > 12 & Tmin > -3 & Tmax > 18 ~ "k",
        Tclim1 %in% c("BW","BS0","BS1") & Tmean > 12 & Tmin > -3  ~ "k'",
        Tclim1 %in% c("BW","BS0","BS1")                             ~ "k''",
        # (A)C
        Tclim1 == "(A)C"                                            ~ "a",
        # C y D
        Tclim1 %in% c("C","D") & Tmean > 12 & Tmax > 22 ~ "a",
        Tclim1 %in% c("C","D") & Tmean > 12              ~ "b",
        Tclim1 %in% c("C","D") & md > 3                  ~ "b'",
        Tclim1 %in% c("C","D") & Tmin > -38              ~ "c",
        Tclim1 %in% c("C","D")                            ~ "d",
        TRUE ~ NA_character_
      )
    )
  
  # ── 10. Cuarta letra — Tclim4 (g: máximo térmico anticipado) ──────────────
  # g = el mes más cálido cae antes del solsticio de verano
  # Necesitamos temperatura mensual para comparar meses específicos
  tmax_mes <- temp_mensual |>
    group_by(Station) |>
    summarise(
      mes_mas_calido = Month[which.max(Tc_m)],
      Tc_max5 = max(Tc_m[Month %in% 1:5]),
      Tc_max6_10 = max(Tc_m[Month %in% 6:10]),
      Tc_max7_11 = max(Tc_m[Month %in% 7:11]),
      Tc_max12_4 = max(Tc_m[Month %in% c(12, 1:4)]),
      .groups = "drop"
    )
  
  base <- base |>
    left_join(tmax_mes, by = "Station") |>
    mutate(
      Tclim4 = case_when(
        Hemi == "N" & Tc_max5 > Tc_max6_10  ~ "g",
        Hemi == "S" & Tc_max7_11 > Tc_max12_4 ~ "g",
        TRUE ~ NA_character_
      )
    )
  
  # ── 11. Quinta letra — Tclim5 (w'': canícula) ─────────────────────────────
  base <- base |>
    mutate(
      Tclim5 = case_when(
        Tclim1 %in% c("A","A(C)") & Pver > Ptot - Pver & Hemi == "N" &
          ((Pmm7 < Pmm6 & Pmm7 < Pmm8) | (Pmm8 < Pmm7 & Pmm8 < Pmm6)) ~ "w''",
        Tclim1 %in% c("A","A(C)") & Pver > Ptot - Pver & Hemi == "S" &
          ((Pmm1 < Pmm12 & Pmm1 < Pmm2) | (Pmm2 < Pmm1 & Pmm2 < Pmm3)) ~ "w''",
        TRUE ~ NA_character_
      )
    )
  
  # ── 12. Corrección E polar en Tclim1 y Tclim3=H ──────────────────────────
  base <- base |>
    mutate(
      Tclim3 = if_else(Tclim1 == "E" & Osc < 10, "H", Tclim3),
      Tclim1 = case_when(
        Tclim1 == "E" & Tmin >  0 & Tmean > -2 ~ "E(T)C",
        Tclim1 == "E" & Tmin <= 0 & Tmean > -2 ~ "E(T)",
        Tclim1 == "E" & Tmax <  0               ~ "EF",
        Tclim1 == "E" & Tmin > -10              ~ "EM",
        TRUE ~ Tclim1
      )
    )
  
  # ── 13. Clave climática final ─────────────────────────────────────────────
  base <- base |>
    mutate(
      Clim = paste0(
        Tclim1,
        if_else(is.na(Tclim2), "", Tclim2),
        if_else(is.na(Tclim3), "", Tclim3),
        Ampli,
        if_else(is.na(Tclim4), "", Tclim4),
        if_else(is.na(Tclim5), "", Tclim5)
      )
    )
  
  # ── 14. Selección y renombre final de columnas ────────────────────────────
  sal0 <- base |>
    select(
      Station,
      Pmm   = Ptot,
      Pmin,
      Pmax,
      Tmean,
      Tmin,
      Tmax,
      Pver,
      Regi,
      P_T,
      Osc,
      PPInv,
      Hum   = graph,
      Hemi,
      Tclim1,
      Tclim2,
      Tclim3,
      Tclim4,
      Tclim5,
      Ampli,
      Clim
    )
  
  return(sal0)
}