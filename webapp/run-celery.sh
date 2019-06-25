#!/bin/bash
celery worker -A app.celery --loglevel=info
