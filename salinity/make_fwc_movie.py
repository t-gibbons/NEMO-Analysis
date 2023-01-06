import cv2
import glob

image_folder = '/project/6007519/weissgib/plotting/figs/freshwater_content/'

runid = ['EPM015', 'EPM101', 'EPM102']

for r in runid:
    img_format = r+'_fwc_200m_*.png'

    files = glob.glob(image_folder+img_format)
    files.sort()

    vid = r+'fwc_monthly.avi'

    frame = cv2.imread(files[0])

    height, width, layers = frame.shape

    video = cv2.VideoWriter(vid, 0, 10, (width, height))

    for image in files:
        video.write(cv2.imread(image))

    cv2.destroyAllWindows()
    video.release()
