import cv2
import numpy as np

# P1S1
img = cv2.imread('../Figure6/P1S1_raw.jpg')
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
ret, thresh = cv2.threshold(gray,20, 255, 0)
contours_punch, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_TC89_L1)
cv2.drawContours(img, contours_punch, -1, (0, 255, 255), 3)
cv2.namedWindow("Resized_Window", cv2.WINDOW_NORMAL)
cv2.resizeWindow("Resized_Window", 800, 600)
cv2.imshow('Resized_Window', img)
cv2.waitKey(0)
cv2.destroyAllWindows()

img = cv2.imread('../Figure6/P1S1_raw.jpg')
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
ret, thresh = cv2.threshold(gray,158.5, 255,cv2.THRESH_BINARY_INV)
contours_slide, hierarchy = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
cv2.drawContours(img, contours_slide, -1, (0, 255, 255), 2)
cv2.namedWindow("Resized_Window", cv2.WINDOW_NORMAL)
cv2.resizeWindow("Resized_Window", 1400, 500)
cv2.imshow('Resized_Window', img)
cv2.waitKey(0)
cv2.destroyAllWindows()

fw = open('../Figure6/P1S1_raw_contours_punch.csv',"w")
for i in range(len(contours_punch)):
    contour_sel = contours_punch[i]
    for j in range(len(contour_sel)):
        fw.write(str(contour_sel[j][0][0])+","+str(contour_sel[j][0][1])+","+str(i)+"\n")
fw.close()

fw = open('../Figure6/P1S1_raw_contours_slide.csv',"w")
for i in range(len(contours_slide)):
    contour_sel = contours_slide[i]
    for j in range(len(contour_sel)):
        fw.write(str(contour_sel[j][0][0])+","+str(contour_sel[j][0][1])+","+str(i)+"\n")
fw.close()

#P1S2
img = cv2.imread('../Figure6/P1S2_raw.jpg')
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
ret, thresh = cv2.threshold(gray,120, 255, 0)
contours_punch, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_TC89_L1)
cv2.drawContours(img, contours_punch, -1, (0, 255, 255), 3)
cv2.imshow('image', img)
cv2.waitKey(0)
cv2.destroyAllWindows()

img = cv2.imread('../Figure6/P1S2_raw.jpg')
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
ret, thresh = cv2.threshold(gray,161.5, 255,cv2.THRESH_BINARY_INV)
contours_slide, hierarchy = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
cv2.drawContours(img, contours_slide, -1, (0, 255, 255), 2)
cv2.imshow('image', img)
cv2.waitKey(0)
cv2.destroyAllWindows()

fw = open('../Figure6/P1S2_raw_contours_punch.csv',"w")
for i in range(len(contours_punch)):
    contour_sel = contours_punch[i]
    for j in range(len(contour_sel)):
        fw.write(str(contour_sel[j][0][0])+","+str(contour_sel[j][0][1])+","+str(i)+"\n")
fw.close()

fw = open('../Figure6/P1S2_raw_contours_slide.csv',"w")
for i in range(len(contours_slide)):
    contour_sel = contours_slide[i]
    for j in range(len(contour_sel)):
        fw.write(str(contour_sel[j][0][0])+","+str(contour_sel[j][0][1])+","+str(i)+"\n")
fw.close()

#P1S3
img = cv2.imread('../Figure6/P1S3_raw.jpg')
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
ret, thresh = cv2.threshold(gray,50, 255, 0)
contours_punch, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_TC89_L1)
cv2.drawContours(img, contours_punch, -1, (0, 255, 255), 3)
cv2.imshow('image', img)
cv2.waitKey(0)
cv2.destroyAllWindows()

img = cv2.imread('../Figure6/P1S3_raw.jpg')
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
ret, thresh = cv2.threshold(gray,155, 255,cv2.THRESH_BINARY_INV)
contours_slide, hierarchy = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
cv2.drawContours(img, contours_slide, -1, (0, 255, 255), 2)
cv2.imshow('image', img)
cv2.waitKey(0)
cv2.destroyAllWindows()

fw = open('../Figure6/P1S3_raw_contours_punch.csv',"w")
for i in range(len(contours_punch)):
    contour_sel = contours_punch[i]
    for j in range(len(contour_sel)):
        fw.write(str(contour_sel[j][0][0])+","+str(contour_sel[j][0][1])+","+str(i)+"\n")
fw.close()

fw = open('../Figure6/P1S3_raw_contours_slide.csv',"w")
for i in range(len(contours_slide)):
    contour_sel = contours_slide[i]
    for j in range(len(contour_sel)):
        fw.write(str(contour_sel[j][0][0])+","+str(contour_sel[j][0][1])+","+str(i)+"\n")
fw.close()

#P2S1
import cv2
import numpy as np
img = cv2.imread('../Figure6/P2S1_raw.jpg')
gray = cv2.bitwise_not(cv2.cvtColor(img, cv2.COLOR_BGR2GRAY))
# cv2.namedWindow("Resized_Window", cv2.WINDOW_NORMAL)
# cv2.resizeWindow("Resized_Window", 800, 600)
# cv2.imshow('Resized_Window', gray)
# cv2.waitKey(0)

# img_gray3 = np.zeros((2939,3609))
# for i in range(0,2939):
#     for j in range(0,3609):
#         img_gray3[i][j] = np.max((img[:,:,0][i][j], img[:,:,1][i][j], img[:,:,2][i][j]))
#      print(i)
# img_gray3 = img_gray3.astype(np.uint8)
# cv2.namedWindow("Resized_Window", cv2.WINDOW_NORMAL)
# cv2.resizeWindow("Resized_Window", 800, 600)
# cv2.imshow('Resized_Window', img_gray3)
# cv2.waitKey(0)


ret, thresh = cv2.threshold(gray,80.5, 255,cv2.THRESH_BINARY)
contours_punch, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE,cv2.CHAIN_APPROX_TC89_L1)
cv2.drawContours(img, contours_punch, -1, (0, 255, 255), 3)
cv2.namedWindow("Resized_Window", cv2.WINDOW_NORMAL)
cv2.resizeWindow("Resized_Window", 800, 600)
cv2.imshow('Resized_Window', img)
cv2.waitKey(0)
cv2.destroyAllWindows()

img = cv2.imread('../Figure6/P2S1_raw.jpg')
gray = (cv2.cvtColor(img, cv2.COLOR_BGR2GRAY))
ret, thresh = cv2.threshold(gray,197, 255,cv2.THRESH_BINARY)
contours_slide, hierarchy = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
cv2.drawContours(img, contours_slide, -1, (0, 255, 255), 2)
cv2.namedWindow("Resized_Window", cv2.WINDOW_NORMAL)
cv2.resizeWindow("Resized_Window", 800, 600)
cv2.imshow('Resized_Window', img)
cv2.waitKey(0)
cv2.destroyAllWindows()

#circled
img = cv2.imread('../Figure6/P2S1_circle.jpg')
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
ret, thresh = cv2.threshold(gray,40, 255, 0)
contours_punch, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_TC89_L1)
cv2.drawContours(img, contours_punch, -1, (0, 255, 255), 3)
cv2.namedWindow("Resized_Window", cv2.WINDOW_NORMAL)
cv2.resizeWindow("Resized_Window", 800, 600)
cv2.imshow('Resized_Window', img)
cv2.waitKey(0)
cv2.destroyAllWindows()

img = cv2.imread('../Figure6/P2S1_circle.jpg')
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
ret, thresh = cv2.threshold(gray,167, 255,cv2.THRESH_BINARY_INV)
contours_slide, hierarchy = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
cv2.drawContours(img, contours_slide, -1, (0, 255, 255), 2)
cv2.namedWindow("Resized_Window", cv2.WINDOW_NORMAL)
cv2.resizeWindow("Resized_Window", 800, 600)
cv2.imshow('Resized_Window', img)
cv2.waitKey(0)
cv2.destroyAllWindows()

fw = open('../Figure6/P2S1_raw_contours_punch.csv',"w")
for i in range(len(contours_punch)):
    contour_sel = contours_punch[i]
    for j in range(len(contour_sel)):
        fw.write(str(contour_sel[j][0][0])+","+str(contour_sel[j][0][1])+","+str(i)+"\n")
fw.close()

fw = open('../Figure6/P2S1_raw_contours_slide.csv',"w")
for i in range(len(contours_slide)):
    contour_sel = contours_slide[i]
    for j in range(len(contour_sel)):
        fw.write(str(contour_sel[j][0][0])+","+str(contour_sel[j][0][1])+","+str(i)+"\n")
fw.close()

#P3S1
#circled
img = cv2.imread('../Figure6/P3S1_circle.jpg')
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
ret, thresh = cv2.threshold(gray,40, 255, 0)
contours_punch, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_TC89_L1)
cv2.drawContours(img, contours_punch, -1, (0, 255, 255), 3)
cv2.namedWindow("Resized_Window", cv2.WINDOW_NORMAL)
cv2.resizeWindow("Resized_Window", 800, 600)
cv2.imshow('Resized_Window', img)
cv2.waitKey(0)
cv2.destroyAllWindows()

img = cv2.imread('../Figure6/P3S1_circle.jpg')
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
ret, thresh = cv2.threshold(gray,187.8, 255,cv2.THRESH_BINARY_INV)
contours_slide, hierarchy = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
cv2.drawContours(img, contours_slide, -1, (0, 255, 255), 2)
cv2.namedWindow("Resized_Window", cv2.WINDOW_NORMAL)
cv2.resizeWindow("Resized_Window", 800, 600)
cv2.imshow('Resized_Window', img)
cv2.waitKey(0)
cv2.destroyAllWindows()

fw = open('../Figure6/P3S1_raw_contours_punch.csv',"w")
for i in range(len(contours_punch)):
    contour_sel = contours_punch[i]
    for j in range(len(contour_sel)):
        fw.write(str(contour_sel[j][0][0])+","+str(contour_sel[j][0][1])+","+str(i)+"\n")
fw.close()

fw = open('../Figure6/P3S1_raw_contours_slide.csv',"w")
for i in range(len(contours_slide)):
    contour_sel = contours_slide[i]
    for j in range(len(contour_sel)):
        fw.write(str(contour_sel[j][0][0])+","+str(contour_sel[j][0][1])+","+str(i)+"\n")
fw.close()

import cv2
import numpy as np

# 读取原始图像
img = cv2.imread('../Figure6/P1S2_raw.jpg')
w,h,c = img.shape
# 定义原始图像中的四个点和目标矩形的四个点
src_pts = np.float32([[141, 131], [480, 159], [493, 630], [64, 601]])
dst_pts = np.float32([[0, 0], [280, 0], [280, 479], [0, 479]])
# 计算透视变换矩阵
M = cv2.getPerspectiveTransform(src_pts, dst_pts)
# 应用透视变换
warped_img = cv2.warpPerspective(img, M,(w,h))
# 显示结果
# cv2.imshow("Original Image", img)
cv2.namedWindow("Resized_Window", cv2.WINDOW_NORMAL)
cv2.resizeWindow("Resized_Window", 800, 600)
cv2.imshow('Resized_Window', warped_img)
cv2.waitKey(0)
cv2.destroyAllWindows()

#HE images
import cv2
import numpy as np

img = cv2.imread('C:/Users/zhang/Downloads/2023_04_17__high_quality_HE_images/2023_04_17__0033/P1.tif')
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
ret, thresh = cv2.threshold(gray,23.8, 255,cv2.THRESH_BINARY_INV)
contours_slide, hierarchy = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
# contours_punch, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_TC89_L1)
cv2.drawContours(img, contours_slide, -1, (0, 255, 255), 2)
cv2.namedWindow("Resized_Window", cv2.WINDOW_NORMAL)
cv2.resizeWindow("Resized_Window", 1600, 800)
cv2.imshow('Resized_Window', img)
cv2.waitKey(0)
cv2.destroyAllWindows()

fw = open('../Figure6/P1_contours_HE.csv',"w")
for i in range(len(contours_slide)):
    contour_sel = contours_slide[i]
    for j in range(len(contour_sel)):
        fw.write(str(contour_sel[j][0][0])+","+str(contour_sel[j][0][1])+","+str(i)+"\n")
fw.close()


img = cv2.imread('C:/Users/zhang/Downloads/2023_04_17__high_quality_HE_images/2023_04_17__0033/P1.png')
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
# ret, thresh = cv2.threshold(gray,10.8, 255,cv2.THRESH_BINARY_INV)
# contours_slide, hierarchy = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
ret, thresh = cv2.threshold(gray,220, 255, 0)
contours_punch, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_TC89_L1)
cv2.drawContours(img, contours_punch, -1, (0, 255, 255), 3)
# cv2.drawContours(img, contours_slide, -1, (0, 255, 255), 2)
cv2.namedWindow("Resized_Window", cv2.WINDOW_NORMAL)
cv2.resizeWindow("Resized_Window", 1600, 800)
cv2.imshow('Resized_Window', img)
cv2.waitKey(0)
cv2.destroyAllWindows()



img = cv2.imread('../Figure6/P1S1_OutCountour_line-01.jpg')
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
ret, thresh = cv2.threshold(gray,191.5, 255,cv2.THRESH_BINARY_INV)
contours_slide, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
# ret, thresh = cv2.threshold(gray,220, 255, 0)
# contours_punch, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_TC89_L1)
cv2.drawContours(img, contours_slide, -1, (0, 0, 255), 3)
# cv2.drawContours(img, contours_slide, -1, (0, 255, 255), 2)
cv2.namedWindow("Resized_Window", cv2.WINDOW_NORMAL)
cv2.resizeWindow("Resized_Window", 1600, 800)
cv2.imshow('Resized_Window', img)
cv2.waitKey(0)
cv2.destroyAllWindows()

fw = open('../Figure6/P1S1_contours_region.csv',"w")
for i in range(len(contours_slide)):
    contour_sel = contours_slide[i]
    for j in range(len(contour_sel)):
        fw.write(str(contour_sel[j][0][0])+","+str(contour_sel[j][0][1])+","+str(i)+"\n")
fw.close()

#liwei
img = cv2.imread('C:/Users/zhang/Downloads/LiWei20230430004015.jpg')
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
# ret, thresh = cv2.threshold(gray,40.5, 255,cv2.THRESH_BINARY_INV)
# contours_slide, hierarchy = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
ret, thresh = cv2.threshold(gray,180, 255, 0)
contours_punch, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_TC89_L1)
cv2.drawContours(img, contours_punch, -1, (0, 120, 255), 3)
# cv2.drawContours(img, contours_slide, -1, (0, 255, 255), 2)
cv2.namedWindow("Resized_Window", cv2.WINDOW_NORMAL)
cv2.resizeWindow("Resized_Window", 1000, 600)
cv2.imshow('Resized_Window', img)
cv2.waitKey(0)
cv2.destroyAllWindows()


import os
os.chdir("C:/wps/0ther/DongZhen")

from bs4 import BeautifulSoup
with open('054200000010noCapID.xml',encoding='UTF-8') as fp:
    soup = BeautifulSoup(fp, 'xml')

def tag_contains_keyword(tag, keyword):
    return tag.name and keyword in tag.name

keyword = 'Shape_'
matching_tags = soup.find_all(lambda tag: tag_contains_keyword(tag, keyword))

fw = open("054200000010noCapID_loc.txt","w",encoding='UTF-8')
i=1
for tag in matching_tags:
    X_all = [z.text for z in tag.find_all(lambda tag: tag_contains_keyword(tag, "X_"))]
    Y_all = [z.text for z in tag.find_all(lambda tag: tag_contains_keyword(tag, "Y_"))]
    fw.write(str(i)+'\t'+";".join(X_all)+'\t'+";".join(Y_all)+"\n")
    i=i+1
fw.close()
matching_tags_Xcali = soup.find_all(lambda tag: tag_contains_keyword(tag, "X_CalibrationPoint"))
matching_tags_Ycali = soup.find_all(lambda tag: tag_contains_keyword(tag, "Y_CalibrationPoint"))

fw = open("LCM/054200000010noCapID_cali.txt","w",encoding='UTF-8')
for i in range(len(matching_tags_Xcali)):
    fw.write(matching_tags_Xcali[i].text+"\t"+matching_tags_Ycali[i].text+'\n')
fw.close()
