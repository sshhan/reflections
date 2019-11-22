#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'shangshuhan'


from flask_cors import CORS
import jieba
import jieba.analyse
from flask import Flask, request, jsonify, render_template
import jieba.posseg as pseg
from face import allowed_file, init_result, FaceCompare, FaceCharacter

app = Flask(__name__)
CORS(app, resources=r'/*')

@app.route('/')
@app.route('/home')
def home():
    return render_template('index.html')

@app.route('/preproc', methods=["GET"])
def pre_process():
    return render_template('preproc.html')

@app.route('/keyExt', methods=["GET"])
def keyword_extraction():
    return render_template('keyExt.html')

@app.route('/summary', methods=["GET"])
def summary():
    return render_template('text_Summarization.html')


@app.route('/installation', methods=["GET"])
def installation():
    return render_template('installation.html')

@app.route('/others', methods=["GET"])
def others():
    return render_template('others.html')


@app.route("/ai/nlp/ChineseTokenize", methods=['POST'])
def chinese_tokenize():
    # 中文分词
    text = request.form['text']
    seg_list = jieba.cut(text, cut_all=False)
    word = "/".join(seg_list)
    result = {
        "result": word
    }
    result = {str(key): value for key, value in result.items()}
    return jsonify(result=result)

@app.route("/ai/nlp/ExtratTags", methods=['POST'])
def extract_tags():
    # 关键词提取
    text = request.form['text']
    topK = request.form['topK']
    word = jieba.analyse.extract_tags(text, topK=int(topK))
    result = {
        "result": word
    }
    result = {str(key): value for key, value in result.items()}

    return jsonify(result=result)

@app.route("/ai/nlp/PosTags", methods=['POST'])
def pos_tags():
    # 词性标注
    text = request.form['text']
    word = []
    for w, tag in pseg.lcut(text):
        word.append('%s(%s)' % (w, tag))
    result = {
        "result": word
    }
    result = {str(key): value for key, value in result.items()}
    return jsonify(result=result)

@app.route("/ai/cv/FaceCompare", methods=["POST", "GET"])
def cvFaceCompare():
    result = {}
    result = init_result(result)
    if request.method == "POST":
        result["methods"] = "POST"

        file1 = request.files['file1']
        file2 = request.files['file2']

        if file1 and allowed_file(file1.filename) and file2 and allowed_file(file2.filename):
            result["success"] = True
            FaceCompare(result, file1, file2)

        else:
            result['code'] = 400
            result['data']['message'] = 'There is no file or file extension is not allowed'
    else:
        result['code'] = 404
        result['data']['message'] = 'please use post method'

    return jsonify(result)

@app.route("/ai/cv/FaceCharacter", methods=["POST", "GET"])
def cvFaceCharacter():
    result = {}
    result = init_result(result)
    if request.method == "POST":
        result["methods"] = "POST"

        file1 = request.files['file']

        if file1 and allowed_file(file1.filename):
            result["success"] = True
            FaceCharacter(result, file1)
        else:
            result['code'] = 400
            result['data']['message'] = 'There is no file or file extension is not allowed'
    else:
        result['code'] = 404
        result['data']['message'] = 'please use post method'

    return jsonify(result)

if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5000, debug=True)

