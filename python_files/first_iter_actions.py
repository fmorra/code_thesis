def first_iter_actions(base_path, base_e_path, new_e_iter_path):
    import numpy as np
    import os

    no_iter_path = os.path.join(base_path, 'no_iter_e.dat')
    if not os.path.exists(no_iter_path):
        with open(base_e_path, "r") as f:
            data = f.readlines()
            data = [line.strip() for line in data]
            data = [line.strip().replace("  ", " ").replace(" ", ",") for line in data]
            data = [np.array([eval(i) for i in line.split(",")[:]]) for line in data]
            data = np.array(data)
            new_column = np.zeros((len(data), 1))
            data = np.column_stack((data, new_column))
        data = np.matrix(data)
        # with open(new_e_iter_path, 'wb') as f_save:
        #     np.savetxt(f_save, data, fmt="%s")
        with open(no_iter_path, 'wb') as f_save:
            np.savetxt(f_save, data, fmt="%s")
    else:
        os.remove(base_e_path)
        with open(no_iter_path, "r") as f:
            data = f.readlines()
            data = [line.strip() for line in data]
        data = np.matrix(data)
        # with open(new_e_iter_path, 'wb') as f_save:
        #     np.savetxt(f_save, data, fmt="%s")
        with open(base_e_path, 'wb') as f_save:
            np.savetxt(f_save, data, fmt="%s")

    return no_iter_path
