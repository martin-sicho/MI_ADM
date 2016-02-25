from django.http import HttpResponse
from django.shortcuts import render

def home(req):
    return HttpResponse('Running...')
