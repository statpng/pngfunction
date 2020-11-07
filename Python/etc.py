def png_unzip(zip_path):
    from zipfile import ZipFile
    import os
    import re
    # zip_path = "./porto-seguro-safe-driver-prediction/porto-seguro-safe-driver-prediction.zip"
    zf = ZipFile(zip_path, 'r')
    dir_path = re.sub(".zip", "", zip_path)
    if not os.path.isdir(dir_path):
        os.mkdir(dir_path)
    zf.extractall(dir_path)
    zf.close()
