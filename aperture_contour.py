#coding:utf-8
import numpy as np

#4方向すべて試してみて進める方向があるか確認
def go_point(point, points):
    directions = [(0.5, 0), (0, 0.5), (-0.5, 0), (0, -0.5)]
    #各方向に進めるか確認
    for direction in directions:
        point_cand = point + direction
        is_point = np.all(points == point_cand, axis=1)
        #もし進めるなら
        if np.sum(is_point) != 0:
            #もとの点を削除
            indexes = np.where(is_point)
            points = np.delete(points, indexes[0][0], 0)
            return point_cand, points
    #どこにも進めない場合、Noneを返す
    else:
        return None, points

def make_contours(aperture):
    #apertureの各点の座標を求める
    h, w = np.where(aperture == 1)
    pixs = np.vstack((w, h)).T
    #各apertureの端の点と辺の点を求める
    TR = pixs + np.array([-0.5, -0.5])
    TL = pixs + np.array([0.5, -0.5])
    BR = pixs + np.array([-0.5, 0.5])
    BL = pixs + np.array([0.5, 0.5])
    R = pixs + np.array([-0.5, 0])
    L = pixs + np.array([0.5, 0])
    T = pixs + np.array([0, -0.5])
    B = pixs + np.array([0, 0.5])
    points_edge = np.vstack((TR, TL, BR, BL))
    points_side = np.vstack((T, B, R, L))
    #重複を削除したものを作成
    points_edge_unique = np.unique(points_edge, axis=0)
    points_side_unique = np.unique(points_side, axis=0)
    #apertureを構成する点ではないものを削除
    points_edge_all = np.array([point for point in points_edge_unique if np.sum(np.all(points_edge == point, axis=1)) != 4])
    points_side_all = np.array([point for point in points_side_unique if np.sum(np.all(points_side == point, axis=1)) != 2])
    #crossする点は重複カウントする
    #crossする点は2つのpixelのedgeを共有
    points_cross_cand = [point for point in points_edge_all if np.sum(np.all(points_edge == point, axis=1)) == 2]
    #crossする点は近傍4辺を持つ
    points_cross_all = np.array([point for point in points_cross_cand if np.sum(np.all(np.abs(points_side_all - point) <= 0.5, axis=1)) == 4])
    #まとめる
    if points_cross_all.shape[0] == 0:
        points_all = np.vstack((points_edge_all, points_side_all))
    else:
        points_all = np.vstack((points_edge_all, points_side_all, points_cross_all))
    #lineを作る点を並べる
    lines = []
    #apertureは一つとは限らないので、全てのaperutreを経由できているか確認
    while points_all.shape[0] != 0:
        point = points_all[0]
        contour = []
        #一つのaperuter内の輪郭を引く
        while point is not None:
            contour.append(point)
            point, points_all = go_point(point, points_all)
        #折れ曲がり点のみ抽出
        contour_arr = np.array(contour)
        diffs = np.diff(contour_arr, 2, axis=0)
        diffs = np.vstack((np.zeros(2), diffs, np.zeros(2)))
        contour_point = np.vstack((contour[0], contour_arr[np.all(diffs != 0, axis=1)], contour[0]))
        lines.append(contour_point)
    return lines

def draw_contours(canvas, aperture, **kwargs):
    points_list = make_contours(aperture)
    for points in points_list:
        for i in range(len(points) - 1):
            canvas.plot((points[i][0], points[i+1][0]), (points[i][1], points[i+1][1]),  "-", **kwargs)
    return canvas
