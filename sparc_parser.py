def get_galaxy_geometry_from_table(galaxy_name):
    # Очистка имени и чтение Table1.mrt (из sparc_to_correct.py)
    search_name = galaxy_name.replace(' ', '')
    with open("sparc_to_correct.py", "r", encoding="utf-8") as f:
        lines = f.readlines()
        
    for line in lines:
        parts = line.strip().split()
        if not parts or parts[0] != search_name: continue
            
        # Логика из sparc_to_correct.py: 
        # r_disk (col 12), тип (col 2), расчет z_c и m_bh
        r_disk = float(parts[11])
        g_type = int(parts[1])
        z_c = r_disk * (0.25 if g_type >= 9 else 0.08 if g_type <= 3 else 0.12)
        m_bh = 0.001 if "3198" in search_name else 0.0001
        return r_disk, z_c, m_bh
            
    return 4.0, 0.35, 0.001 # Дефолт
